import astropy.units as u
from astropy.time import Time, TimeDelta

from radio_pulsar_nav.config import fast, kens, parkes
from radio_pulsar_nav.psrcat import calculate_observable_pulsars, load_catalogue
from radio_pulsar_nav.visualise import (
    plot_antenna_coverage,
    plot_flux_density,
    plot_pulsar_position,
)


def test_main():
    antenna = kens
    cat, flux_cat = load_catalogue(antenna.centre_frequency)
    t = Time.now()
    horiz = True
    for integ_time in [
        1e10 * 60,
        1e10 * 60,
    ]:  # 60 * 60, 6*60*60, 24*60*60]:
        fig_dens = plot_flux_density(cat, flux_cat)
        vis = calculate_observable_pulsars(cat, flux_cat, antenna, t, integ_time, above_horiz=horiz)
        horiz = not horiz
        col = ["rgba(238,118,0,0.7)" if psr else "rgba(104,131,139,0.4)" for psr in vis]
        fig_dens.add_hline(
            y=antenna.min_observable_flux_density(integ_time),
            line_dash="dash",
            label=dict(text=antenna.name, textposition="end"),
        )
        fig_dens.update_traces(marker=dict(color=col))
        fig_dens.update_layout(
            title=dict(
                text=f"Pulsar observability over {integ_time}s integration time",
                xanchor="center",
                yanchor="top",
            )
        )
        fig_dens.show()

        fig_pos = plot_pulsar_position(cat)
        fig_pos.update_traces(marker=dict(color=col))
        fig_pos.update_layout(
            title=dict(
                text=f"Pulsar Positioning over {integ_time}s integration time",
                xanchor="center",
                yanchor="top",
            )
        )
        plot_antenna_coverage(fig_pos, antenna, t)
        fig_pos.show()
