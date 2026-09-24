"""Independent-cloud superposition, common origin and simultaneous start."""
from dataclasses import asdict
import numpy as np
from cloud import Gas, Atmosphere, Vessel, Options, PrimaryCloud
from secondary import GasFeed, PlumeOptions, SecondaryCloud


def validate_pair(primary, secondary):
    from flammable_mass import limits, validate as validate_mass
    validate_mass(primary);validate_mass(secondary)
    a_limits,b_limits=limits(primary),limits(secondary)
    if a_limits != b_limits:
        raise ValueError("Для общей массы задайте одинаковые НКПР и ВКПР в обеих вкладках")
    from exposure import defaults, validate
    defaults(primary);defaults(secondary);validate(primary);validate(secondary)
    if primary['toxicity'] != secondary['toxicity']:
        raise ValueError('Для суммы задайте одинаковые PCt50 и LCt50 в обеих вкладках')
    for section, cls, ignored in [('gas', Gas, set()), ('atmosphere', Atmosphere, {'alpha_basis'})]:
        a, b = asdict(cls(**primary[section])), asdict(cls(**secondary[section]))
        different = [k for k in a if k not in ignored and a[k] != b[k]]
        if different:
            raise ValueError('Для суммы нужны одинаковые вещество и погода: '
                             + section + ': ' + ', '.join(different)
                             + '. Используйте кнопку копирования вещества и погоды.')
    if primary.get('section_height_m', 0.) != secondary.get('section_height_m', 0.):
        raise ValueError('Для общей карты задайте одинаковую высоту сечения обоих облаков')


class CombinedCloud:
    def __init__(self, primary, secondary):
        self.primary, self.secondary = primary, secondary
        self.time_end = min(float(primary.solution.t[-1]), secondary.options.observation_s)
        self.times = np.unique(np.append(primary.times[primary.times <= self.time_end], self.time_end))

    @classmethod
    def from_inputs(cls, a, b):
        validate_pair(a, b)
        p = PrimaryCloud(Gas(**a['gas']), Vessel(**a['vessel']), Atmosphere(**a['atmosphere']), Options(**a.get('options', {})))
        s = SecondaryCloud(Gas(**b['gas']), GasFeed(**b['feed']), Atmosphere(**b['atmosphere']), PlumeOptions(**b.get('plume_options', {})))
        p.run(); s.run()
        return cls(p, s)

    def concentration(self, t, x, y=0., z=0.):
        if not np.isfinite(t) or not 0 <= t <= self.time_end:
            raise ValueError('Время вне общего рассчитанного интервала')
        a = self.primary.concentration(t, x, y, z)
        b = self.secondary.concentration(x, y, z, t)
        return a, b, a + b

    def maximum(self, x, y=0., z=0., times=None):
        """Sample maximum of the sum; never sum separate maxima."""
        times = self.times if times is None else times
        peak = self.concentration(0., x, y, z)[2]
        for t in times:
            peak = np.maximum(peak, self.concentration(float(t), x, y, z)[2])
        return peak

    def receptor_times(self, x, extra=()):
        times = list(self.times) + [t for t in extra if 0 <= t <= self.time_end]
        if 0 <= x <= self.secondary.x[-1]:
            arrival = float(self.secondary.state(x, self.secondary.solution.sol(x))['travel_time_s'])
            for event in (arrival, arrival + self.secondary.feed.duration_s):
                # Resolve both sides of the discontinuous finite-feed gate.
                times.extend(t for t in (np.nextafter(event, -np.inf), event)
                             if 0 <= t <= self.time_end)
        return np.unique(times)
