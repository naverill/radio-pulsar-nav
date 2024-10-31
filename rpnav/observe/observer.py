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
        super().__init__(name, aliases, timescale)
        self.update(
            location=location,
            origin=origin,
            time=time,
            itoa_code=itoa_code,
        )

    def to_json(self):
        t = self._time if self._time is not None else Time.now()
        itrf = self._location.get_itrs(t)
        out = {
            self.name: {
                "itrf_xyz": [itrf.x.to_value(u.m), itrf.y.to_value(u.m), itrf.z.to_value(u.m)],
            }
        }
        if self._origin:
            out[self.name]["origin"] = self._origin
        if self._itoa_code:
            out[self.name]["itoa_code"] = self._itoa_code
        return out

    def update(
        self,
        location: EarthLocation = None,
        origin: str = None,
        time: Time = None,
        itoa_code: str = None,
    ):
        if location is not None:
            self._location = location
        if origin is not None:
            self._origin = origin
        if time is not None:
            self._time = time
        if itoa_code is not None:
            self._itoa_code = itoa_code
        try:
            load_observatories(io.StringIO(json.dumps(self.to_json())), overwrite=True)
        except ValueError:
            load_observatories(io.StringIO(json.dumps(self.to_json())))
        return self