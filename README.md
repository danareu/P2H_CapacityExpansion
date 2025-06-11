# P2H_CapacityExpansion
Capacity Expansion Model for H2 and Power Investment Decisions 

P2H_CapacityExpansion is a julia implementation of an input-data scaling capacity expansion modeling framework.

The primary purpose of this model is to serve as a baseline model for a surrogate model. 
Specifically, the following workflow is defined:

![Picture1](https://github.com/user-attachments/assets/20ca0c3c-0c3a-41fc-8e58-17808023e0aa)<br/>


## The Model
The model can be seen as a simplified version of the energy system model [GENeSYS-MOD](https://github.com/GENeSYS-MOD/GENeSYS_MOD.jl). To better make use of Julia's capabilities, the general structure and modeling are inspired by the power sector model [CapacityExpansion.jl](https://github.com/YoungFaithful/CapacityExpansion.jl). <br/>

The dispatch module, which will be approximated by a data-driven approach, can be described as follows:

![dispatch](https://github.com/user-attachments/assets/9713454e-22a6-4bca-a21a-d20a9119c41c)


## The Data
The technology data is based on the EU-funded [Man0EUvRE project](https://man0euvre.eu/). The corresonding data and data sources are available on [git](https://github.com/GENeSYS-MOD/GENeSYS_MOD.data/tree/development/man0evure_scenario_refresh). The annual hydrogen demand for each country was obtained from the European Hydrogen Backbone report [1]. The inflow profiles and hydropower storage capacities were taken from the latest [TYNDP report](https://2024.entsos-tyndp-scenarios.eu/download/). The weekly inflow profiles are based on the weather year 2016. The hourly capacity factors for solar PV, offshore and onshore wind are obtained from renewables.ninja [2,3,4], a widely recognized open-access platform that derives its estimates from NASA’s MERRA-2 reanalysis data. Furthermore, we use hourly electricity demand data using the ENTSO-E API [4]. These time-series are based on the weather year 2020. <br/>



[1] European Hydrogen Backbone. (2022). A European hydrogen infrastructure vision covering 28 countries. https://ehb.eu/files/downloads/ehb-report-220428-17h00-interactive-1.pdf<br/>
[2] I. Staffell, S. Pfenninger, Using bias-corrected reanalysis to simulate current and future wind power output, Energy 114 (2016) 1224–1239. doi:https://doi.org/10.1016/j.energy.2016.08.068.<br/>
[3] Renewables.ninja, https://renewables.ninja, accessed: Jul. 11, 2024.<br/>
[3] S. Pfenninger, I. Staffell, Long-term patterns of european pv output using 30 years of validated hourly reanalysis and satellite data, Energy 114 (2016) 1251–1265. doi:10.1016/j.energy.2016.08.060.<br/>
[4] EnergieID, entsoe-py, https://github.com/EnergieID/entsoe-py, accessed: 2025-01-09 (2025).
