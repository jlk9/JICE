# Written by Joseph Kump (josek97@utexas.edu)
# Constants used throughout JICE models

# Used in single-column thermodynamics:
const V_a    = 0.025
const K_a    = 0.03
const T_mlt0 = 273.15   # K         melting point of freshwater
const ρ_0    = 917.0    # kg/m^3    density of fresh (pure) ice
const ρ_i    = 917.0    # kg/m^3    density of sea ice
const ρ_w    = 1025.0   # kg/m^3    density of seawater (based on common estimate)
const ρ_s    = 330.0    # kg/m^3    density of snow
const c_0    = 2106.0   # J/kg/K    specific heat of fresh ice
const c_w    = 3900.0   # J kg K    specific heat of seawater (estimate)
const c_ocn  = 4218.0   # J kg K    specific heat of the ocean
const c_h    = 0.006    #           heat transfer coefficient
const K_s    = 0.3      # W/m/deg   thermal conductivity of snow
const L_0    = 334000.0 # J / kg    latent heat of fusion of fresh ice
const κ_i    = 1.4      #           extinction coefficient
const μ      = 0.054    # deg/ppt   liquidus ratio between the freezing temperature and salinity of brine
const S_max  = 3.2      # ppt       maximum salinity of sea ice
const ahmax  = 0.3      # m         thickness above which albedo is constant
const α_o    = 0.06     #           ocean albedo

const α_icev    = 0.78  # visible ice albedo for h > ahmax
const α_icei    = 0.36  # near-ir ice albedo for h > ahmax
const α_snowv   = 0.98  # cold snow albedo, visible
const α_snowi   = 0.70  # cold snow albedo, near IR

const dα_mlt    = -0.075 # albedo change for temp change of 1 degree for ice
const dα_mltv   = -0.1   # albedo change for temp change of 1 degree for snow, visible
const dα_mlti   = -0.15  # albedo change for temp change of 1 degree for snow, infrared

const i0vis     = 0.7   # fraction of penetrating solar radiation at top surface

const puny       = 1.0e-11  # For numerical tests
const Tsf_errmax = 0.01     # For numerical test of convergence
const hfrazilmin = 0.05     # m minimum thickness of new frazil ice


# Used in grid cell thermodynamics:
const sea_sal          = 35.0    # typical sea salt in ppt
const dSin0_frazil     = 3.0     # bulk salinity reduction of newly formed frazil
const ice_ref_salinity = 4.0     # reference salinity of sea ice, ppt


# Used for atmospheric variables:
const emissivity = 0.985                # Emissivity of snow or ice
const sbc        = 5.670374419e-8       # Stefan-Boltzman Constant (FIX)
const q_1        = 1.16378e7            # kg / m^3
const q_2        = 5897.8               # K
const κ          = 0.40                 # von Karman constant
const g          = 9.8                  # gravitational acceleration (positive downward)
const z_ref      = 10.0                 # ice reference height
const U_dmin     = 1.0                  # minimum allowable value for |U_a| (since we don't use high frequency)
const z_deg      = 1.0                  # Assumed level height
const C_to_K     = 273.15               # converts Celcius to Kelvin
const λ          = -2.3025850929940455  # log(z_deg / z_ref)
const L_vap      = 2260.0               # Latent heat of vaporization
const L_ice      = 334.0                # Latent heat of fusion