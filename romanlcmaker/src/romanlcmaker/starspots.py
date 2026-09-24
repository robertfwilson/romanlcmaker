import butterpy as bp

import astropy.units as u
from astropy.modeling.models import BlackBody


class Starspots(object):
    def __init__(self, time_array, spot_params):
        self.time_array = time_array
        self.spot_params = spot_params

    def generate_lightcurve(self):
        # Create a starspot model using butterpy
        model = bp.StarspotModel(self.time_array, self.spot_params)
        lightcurve = model.generate_lightcurve()


        Prot = self.spot_params['Prot']
        inc = self.spot_params['inc']
        activity_level = self.spot_params['activity_level']
        activity_cycle_period = self.spot_params['activity_cycle_period']
        activity_cycle_overlap = self.spot_params['activity_cycle_overlap']


        spot_lifetime = self.spot_params['spot_lifetime']
        rot_shear = self.spot_params['rot_shear']
        
        min_lat = self.spot_params['min_lat']
        max_lat = self.spot_params['max_lat']
        star_teff = self.spot_params['star_teff']
        spot_teff = star_teff - self.spot_params['delta_spot_teff']

            # Emerge active regions
        s = bp.Surface()
        
        
        burn_in_time = 10.*Prot

        scale=1.

        regions = s.emerge_regions( ndays=burn_in_time + np.max(self.time_array), 
            activity_level=activity_level*scale, 
            cycle_period=activity_cycle_period,
            cycle_overlap=activity_cycle_overlap,
                min_lat=min_lat, max_lat=max_lat, )
        
        
        

        time_season=self.time_array + burn_in_time

        lightcurve_other = s.evolve_spots(time=time_season,
                                        inclination=inc,
                                        period=Prot,
                                        shear=rot_shear,  
                                        tau_evol=spot_lifetime,
                                        alpha_med=1e-2,
                                        threshold=0.5, )
        
        

        flux_season =  lightcurve_other.flux #CubicSpline( x=lightcurve_other.time - burn_in_time, y=lightcurve_other.flux, )(time_season)


        spot_diff_f146, spot_diff_f213, spot_diff_f087 = get_multiband_rotation_lightcurves(star_fraction=flux_season, 
                                                                                        star_teff=star_teff, 
                                                                                        spot_teff=spot_teff,
                                                                                        scale=1,  band_waves=[1.464, 2.125, 0.869], )




    return lightcurve




def get_multiband_rotation_lightcurves(star_fraction, star_teff, spot_teff, scale=1., 
                                       band_waves=[1.45, 2.1, 0.87], ): #bands = ['F146','F087','F213']):

    spot_fraction = (1.-star_fraction)
    star_fraction_scaled = 1.-spot_fraction

    spot_bb = BlackBody(temperature = spot_teff, )
    star_bb = BlackBody(temperature = star_teff, )

    spot_fluxes = [spot_fraction * (spot_bb(wave*u.um).value/star_bb(wave*u.um).value) for wave in band_waves]

    diff_lightcurves = [scale*(sp_flux+star_fraction - np.mean(sp_flux+star_fraction) ) for sp_flux in spot_fluxes]
    
    return diff_lightcurves
