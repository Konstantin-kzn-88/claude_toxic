"""Experimental single-phase primary cloud, Rostechnadzor Order 385 (2022).
All quantities in SI. Equation numbers follow the supplied PDF, not its cross-links.
No secondary cloud, phase equilibrium, toxic effect or risk calculation.
"""
from dataclasses import dataclass, asdict, fields
from math import pi, sqrt, gamma, log, atan, isfinite
import numpy as np
from scipy.integrate import solve_ivp

R_GAS = 8.3144  # Appendix 1, J/(mol K)
G = 9.81

@dataclass(frozen=True)
class Gas:
    name: str
    molar_mass_kg_mol: float
    adiabatic_index: float
    minimum_valid_temperature_k: float
    phase_and_ideal_gas_assumptions_confirmed: bool

@dataclass(frozen=True)
class Vessel:
    volume_m3: float
    pressure_abs_pa: float
    temperature_k: float

@dataclass(frozen=True)
class Atmosphere:
    temperature_k: float
    surface_temperature_k: float
    pressure_pa: float
    wind_speed_m_s: float
    wind_height_m: float
    roughness_m: float
    stability: str
    alpha_wind: float
    alpha_basis: str
    alpha_height_limit_m: float = 20.0
    air_molar_mass_kg_mol: float = 0.02896
    air_cp_j_kg_k: float = 1005.0

@dataclass(frozen=True)
class Options:
    duration_s: float = 600.0
    output_step_s: float = 0.5
    max_step_s: float = 0.5
    rtol: float = 1e-7
    atol: float = 1e-9


def _positive(name, value):
    if isinstance(value, bool) or not isinstance(value, (float, int)) or not isfinite(value) or value <= 0:
        raise ValueError(f'{name}: требуется конечное положительное число')


def validate(gas, vessel, met, options):
    for obj in (gas, vessel, met, options):
        for field in fields(obj):
            if field.type is float:
                _positive(field.name, getattr(obj, field.name))
    if not gas.name.strip() or not met.alpha_basis.strip():
        raise ValueError('Требуются название вещества и обоснование alpha_wind')
    if gas.phase_and_ideal_gas_assumptions_confirmed is not True:
        raise ValueError('Нужно подтвердить однофазность и применимость идеального газа')
    if gas.adiabatic_index <= 1:
        raise ValueError('Показатель адиабаты должен быть > 1')
    if met.stability not in 'ABCDEF' or len(met.stability) != 1:
        raise ValueError('Класс устойчивости: A, B, C, D, E или F')
    if vessel.pressure_abs_pa < met.pressure_pa:
        raise ValueError('Разрежение в ёмкости не поддерживается')
    if met.air_cp_j_kg_k <= R_GAS / met.air_molar_mass_kg_mol:
        raise ValueError('Теплоёмкость воздуха даёт Cv <= 0')
    if min(vessel.temperature_k, met.temperature_k, met.surface_temperature_k) < gas.minimum_valid_temperature_k:
        raise ValueError('Температура ниже подтверждённого диапазона однофазной модели')
    if options.duration_s > 600:
        raise ValueError('Версия 0.1: t <= 600 с; переменное время осреднения ещё не реализовано')
    if options.duration_s / options.output_step_s > 100000:
        raise ValueError('Слишком много выходных временных точек (максимум 100000)')


def initial_cloud(gas, vessel, met):
    """Eq. (2), (6), (8). Temperature follows ideal-gas/adiabatic closure."""
    mass = gas.molar_mass_kg_mol * vessel.volume_m3 * vessel.pressure_abs_pa / (R_GAS * vessel.temperature_k)
    density = mass / vessel.volume_m3 * (met.pressure_pa / vessel.pressure_abs_pa) ** (1 / gas.adiabatic_index)
    temperature = vessel.temperature_k * (met.pressure_pa / vessel.pressure_abs_pa) ** ((gas.adiabatic_index - 1) / gas.adiabatic_index)
    radius = (mass / (pi * density)) ** (1 / 3)
    if temperature < gas.minimum_valid_temperature_k:
        raise ValueError('Расширение охлаждает газ ниже подтверждённого диапазона; требуется фазовая модель')
    if mass > 500000:
        raise ValueError('Масса > 500 т: вне области этой реализации (п. 11)')
    return dict(mass_kg=mass, density_kg_m3=density, temperature_k=temperature,
                radius_m=radius, height_m=radius, volume_m3=mass / density)


def friction_velocity(met):
    # Eq. (94)-(96), Table 7-6; alpha is supplied separately with its source.
    coefficients = {'A': (-11.4, .10), 'B': (-26., .17), 'C': (-123., .30),
                    'E': (123., .30), 'F': (26., .17)}
    if met.stability == 'D':
        phi = 0.
    else:
        k_l, p = coefficients[met.stability]
        length = k_l * met.roughness_m ** p
        if length > 0:
            phi = -6.9 * met.wind_height_m / length
        else:
            a = (1 - 22 * met.wind_height_m / length) ** .25
            phi = 2 * log((1 + a) / 2) + log((1 + a*a) / 2) - 2 * atan(a) + pi / 2
    denominator = log((met.wind_height_m + met.roughness_m) / met.roughness_m) - phi
    if denominator <= 0:
        raise ValueError('Неположительный знаменатель профиля ветра')
    return .41 * met.wind_speed_m_s / denominator


class PrimaryCloud:
    def __init__(self, gas, vessel, met, options=Options()):
        validate(gas, vessel, met, options)
        self.gas, self.vessel, self.met, self.options = gas, vessel, met, options
        self.initial = initial_cloud(gas, vessel, met)
        self.q = self.initial['mass_kg']
        self.cv_g = R_GAS / gas.molar_mass_kg_mol / (gas.adiabatic_index - 1)
        self.cp_g = self.cv_g * gas.adiabatic_index
        self.cv_a = met.air_cp_j_kg_k - R_GAS / met.air_molar_mass_kg_mol
        self.rho_a = met.pressure_pa * met.air_molar_mass_kg_mol / (R_GAS * met.temperature_k)
        self.u_star = friction_velocity(met)
        self.beta = 1 + met.alpha_wind
        self.h_factor = gamma(1 / self.beta) / self.beta
        self.delta = dict(A=.22, B=.16, C=.11, D=.08, E=.06, F=.04)[met.stability]
        if self.initial['density_kg_m3'] <= self.rho_a:
            raise ValueError('Первоначальное облако не является тяжёлым газом')
        if self.initial['height_m'] > met.alpha_height_limit_m:
            raise ValueError('Начальная высота превышает диапазон выбранного alpha_wind')

    def state(self, y):
        # State = total mass, effective radius, Sy squared, internal energy, centre x.
        total, radius, sy2, energy, x = y
        air = max(total - self.q, 0.)
        heat_capacity = self.q * self.cv_g + air * self.cv_a
        temperature = energy / heat_capacity
        moles = self.q / self.gas.molar_mass_kg_mol + air / self.met.air_molar_mass_kg_mol
        volume = moles * R_GAS * temperature / self.met.pressure_pa  # (209)-(211)
        density = total / volume
        height = volume / (pi * radius**2)
        sz = height / self.h_factor  # (101)
        sz_wind = max(height, .5) / self.h_factor  # footnote 8: only drift speed!
        speed = gamma((1 + self.met.alpha_wind) / self.beta) / gamma(1 / self.beta) * self.met.wind_speed_m_s * (sz_wind / self.met.wind_height_m)**self.met.alpha_wind
        sy2 = max(sy2, 0.)
        core2 = max(radius**2 - sy2, 0.)
        return dict(total_mass_kg=total, radius_m=radius, sy2_m2=sy2, core_radius_m=sqrt(core2),
                    temperature_k=temperature, density_kg_m3=density, height_m=height,
                    sz_m=sz, speed_m_s=speed, centre_x_m=x,
                    centre_concentration_kg_m3=self.q / volume, area_m2=pi * radius**2)

    def rhs(self, t, y):
        s = self.state(y)
        total, radius, sy2, energy, x = y
        temp, rho, height = s['temperature_k'], s['density_kg_m3'], s['height_m']
        cp = (self.q * self.cp_g + max(total - self.q, 0) * self.met.air_cp_j_kg_k) / total
        dt_surface = self.met.surface_temperature_k - temp
        # Literal (114), including the outer square shown in supplied PDF.
        forced = 1.22 * (self.u_star**2 / self.met.wind_speed_m_s)**2 * rho * cp * dt_surface
        natural = 0. if dt_surface <= 0 else 3.5e-3 * (2 * dt_surface)**(2/3) * self.met.pressure_pa / R_GAS * G**(1/3)
        flux = max(forced, natural) if dt_surface > 0 else forced  # (113)-(116)
        w_star = (G * abs(flux) * height / (rho * temp * cp))**(1/3)
        u_t = sqrt(self.u_star**2 + (.2 * w_star)**2)
        ri = G * (rho - self.rho_a) / self.rho_a * height / u_t**2
        phi = sqrt(1 + .8 * ri) / self.beta if ri > 0 else (1 - .6 * ri)**(-.5) / self.beta
        entrainment = .41 * u_t / phi  # (97)-(98)
        radius_dot = 1.15 * sqrt(max(G * height * (1 - self.rho_a / rho), 0))  # (108)
        mass_dot = pi * radius**2 * self.rho_a * entrainment + 2 * pi * radius * height * self.rho_a * .63 * radius_dot
        # (109) multiplied by 2 Sy: regular at Sy=0, no artificial seed.
        sigma_prime = self.delta * (1 + .5e-4 * max(x, 0)) / (1 + 1e-4 * max(x, 0))**1.5
        sy2_dot = 4 * sqrt(2/pi) * s['speed_m_s'] * (s['core_radius_m'] + .5 * sqrt(pi) * sqrt(max(sy2, 0))) * sigma_prime
        energy_dot = mass_dot * self.cv_a * self.met.temperature_k + pi * radius**2 * flux  # (111)
        return [mass_dot, radius_dot, sy2_dot, energy_dot, s['speed_m_s']]

    def run(self):
        y0 = [self.q, self.initial['radius_m'], 0., self.q * self.cv_g * self.initial['temperature_k'], 0.]
        # Stop instead of silently applying unimplemented transitions/invalid physics.
        def core_vanish(t, y): return y[1]**2 - y[2]
        def buoyancy(t, y): return self.state(y)['density_kg_m3'] - self.rho_a
        def phase(t, y): return self.state(y)['temperature_k'] - self.gas.minimum_valid_temperature_k
        def height_limit(t, y): return self.met.alpha_height_limit_m - self.state(y)['height_m']
        def distance_limit(t, y): return 10000 - y[4]
        events = [core_vanish, buoyancy, phase, height_limit, distance_limit]
        for event in events:
            event.terminal, event.direction = True, -1
        sol = solve_ivp(self.rhs, (0., self.options.duration_s), y0, method='DOP853',
                        dense_output=True, events=events, rtol=self.options.rtol,
                        atol=self.options.atol, max_step=self.options.max_step_s)
        if not sol.success:
            raise RuntimeError(sol.message)
        times = np.arange(0., sol.t[-1], self.options.output_step_s)
        times = np.append(times, sol.t[-1])
        states = [self.state(y) for y in sol.sol(times).T]
        names = ['core_disappeared', 'no_longer_heavy', 'temperature_out_of_range',
                 'alpha_height_exceeded', 'distance_10km']
        stop = next((name for name, ev in zip(names, sol.t_events) if len(ev)), 'requested_duration')
        self.solution, self.times, self.states = sol, times, states
        warnings = [
            'Экспериментальная версия: результаты рассеяния не валидированы для проектного применения.',
            'Фазовое равновесие не рассчитывается; однофазность и идеальность подтверждает пользователь.',
            'Максимумы концентрации определяются на конечной временной сетке и в заданном интервале.',
            'Расчёт ограничен 600 с, тяжёлым облаком с ненулевым ядром и заданным диапазоном alpha.',
        ]
        if min(s['height_m'] for s in states) < self.met.roughness_m:
            warnings.append('H < шероховатости: результаты оценочные, п. 21.')
        if stop != 'requested_duration':
            warnings.append(f'Расчёт досрочно остановлен: {stop}; последующее рассеяние не вычислено.')
        return dict(status='experimental_unvalidated', stop_reason=stop,
                    initial=self.initial, warnings=warnings,
                    time_s=times.tolist(), states=states)

    def concentration(self, t, x, y=0., z=0.):
        """Eq. (104)-(106), kg/m³. x/y/z broadcast; t is scalar."""
        if not hasattr(self, 'solution'):
            raise RuntimeError('Сначала вызовите run()')
        if not isfinite(t) or t < 0 or t > self.solution.t[-1]:
            raise ValueError('Время вне рассчитанного интервала')
        x, y, z = np.broadcast_arrays(np.asarray(x, float), np.asarray(y, float), np.asarray(z, float))
        if not all(np.all(np.isfinite(a)) for a in (x, y, z)) or np.any(z < 0):
            raise ValueError('Координаты должны быть конечны, z >= 0')
        s = self.state(self.solution.sol(t))
        d2 = (x - s['centre_x_m'])**2 + y*y
        extra = np.maximum(d2 - s['core_radius_m']**2, 0)
        horizontal = np.exp(-extra / s['sy2_m2']) if s['sy2_m2'] > 0 else (extra == 0).astype(float)
        return s['centre_concentration_kg_m3'] * horizontal * np.exp(-(z / s['sz_m'])**self.beta)

    def maximum_concentration(self, x, y=0., z=0.):
        """Sampled-in-time maximum, NOT a guaranteed global maximum."""
        maximum = self.concentration(self.times[0], x, y, z)
        for t in self.times[1:]:
            maximum = np.maximum(maximum, self.concentration(t, x, y, z))
        return maximum
