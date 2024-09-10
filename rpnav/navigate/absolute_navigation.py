import astropy.units as u
from astropy.coordinates import EarthLocation
from astropy.time import Time
from pint.toa import TOA


class AbsoluteNavigation:
    def __init__(
        self, 
        time_model: TimingModel,
        t0: Time
    ):
        self.t0 = t0
        self.time_model = time_model



    def update_position(
        self, 
        position: EarthLocation,
        pulsar: Pulsar,
        t: Time,
        t_obs: TOA,
    ):  
        """
         ^n is the unit vector to the pulsar with respect
        """
        # calculate predicted arrival tme (t_obs_pred) 
        # delta_t_obs = t_obs_pred - t_obs
        # predicted arrival time at barycentre (t_bary_pred)
        # observed arrival time at barycentre (t_bary)
        # delta_t_bary = t_bary_pred - t_bary  
        # delta_r = (delta_t_bary - delta_t_obs) * SPEED_OF_LIGHT / pulsar_vec

        # pos_icrs = EarthLocation.get_itrs(t).transform_to(ICRS)
        # pulsar_vec = pulsar - pos_icrs

        # delta_t = pulse_arrival_predicted - arrival_measured
        return
        

    def delay(toa: TOA):
        delay = 0 * u.second
        for dc in self.timing_model.DelayComponent_list:
            for df in dc.delay_funcs_component:
                delay += df(toas, delay)
        return delay
