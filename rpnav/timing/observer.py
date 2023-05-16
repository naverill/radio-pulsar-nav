import json

from astropy.coordinates import ITRS, EarthLocation
from astropy.time import Time
from pint.observatory import Observatory


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
        if self.time is None:
            t = Time.now()

        itrf = self.location.get_itrs(t)
        out = {
            self.name: {
                "itrf_xyz": [itrf.x.value, itrf.x.value, itrf.z.value],
            }
        }
        if self.origin:
            out[self.name]["origin"] = self.origin
        if self.itoa_code:
            out[self.name]["itoa_code"] = self.itoa_code
        return out
