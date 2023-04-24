import astropy.units as u
from astropy.time import Time, TimeDelta
from config import fast, kens, parkes
from psrcat import calculate_visible_pulsars, load_catalogue
from visualise import plot_antenna_coverage, plot_flux_density, plot_pulsar_position

antenna = kens
cat, flux_cat = load_catalogue(antenna.centre_freq)
t = Time.now()
t1 = t + TimeDelta(6.0 * 60 * 60 * u.s)

for integ_time in [60, 5 * 60, 15 * 60]:
    fig_dens = plot_flux_density(cat, flux_cat)
    vis = calculate_visible_pulsars(cat, flux_cat, antenna, t, integ_time)
    col = ["rgba(238, 118,0, 0.7)" if psr else "rgba(104,131,139, 0.4)" for psr in vis]
    fig_dens.add_hline(
        y=kens.min_observable_flux_density(integ_time),
        line_dash="dash",
        label=dict(text=antenna.name, textposition="end"),
    )
    fig_dens.update_traces(marker=dict(color=col))
    fig_dens.show()

    fig_pos = plot_pulsar_position(cat)
    fig_pos.update_traces(marker=dict(color=col))
    plot_antenna_coverage(fig_pos, antenna, t)
    plot_antenna_coverage(fig_pos, antenna, t1)
    fig_pos.show()
