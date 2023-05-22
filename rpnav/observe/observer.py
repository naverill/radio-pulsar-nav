import io
import json

import astropy.units as u
from astropy.coordinates import ITRS, EarthLocation
from astropy.time import Time
from pint.observatory import Observatory
from pint.observatory.topo_obs import load_observatories


class Observer(Observatory):
    def __init__(
        self,
        name: str,
        location: EarthLocation,
        time: Time = None,
        origin: str = None,
        aliases: list[str] = [],
        timescale: str = "TCB",
        itoa_code: str = None,
    ):
        self.location = location
        self.origin = origin
        self.time = time
        self.itoa_code = itoa_code
        super().__init__(name, aliases, timescale)

    def to_json(self):
        t = self.time if self.time is not None else Time.now()
        itrf = self.location.get_itrs(t)
        out = {
            self.name: {
                "itrf_xyz": [itrf.x.to_value(u.m), itrf.y.to_value(u.m), itrf.z.to_value(u.m)],
            }
        }
        if self.origin:
            out[self.name]["origin"] = self.origin
        if self.itoa_code:
            out[self.name]["itoa_code"] = self.itoa_code
        return out

    def update(self):
        try:
            load_observatories(io.StringIO(json.dumps(self.to_json())), overwrite=True)
        except ValueError:
            load_observatories(io.StringIO(json.dumps(self.to_json())))
        return self
