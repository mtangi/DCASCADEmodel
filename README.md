The data are divided into three folders:

- Folder DCASCADE_functions contains all scripts and functions used to run D-CASCADE;
- Folder input_data contains the input data necessary to run the two case studies given as examples, the Bega river network and 3S river network  ;
- Folder plot_function reports the function to plot the output data of D-CASCADE.

The main folder also contains two scripts, which show how to run D-CASCADE for the two case studies presented:

- Script main_script_DCASCADE contains all operations necessary to run D-CASCADE on the Bega river network;
- Script main_script_DCASCADE contains all operations necessary to run D-CASCADE on the 3S river network. This case study includes multiple reservoirs in the system, whose features (stored water and sediment volumes, release, flooded area, and others) are dynamically modelled.

A detailed explanation of the model structure and functioning is documented in Tangi, M., Bizzi, S., Fryirs, K., & Castelletti, A. (2022). A dynamic, network scale sediment (dis)connectivity model to reconstruct historical sediment transfer and river reach sediment budgets. Water Resources Research, 58, e2021WR030784. https://doi.org/10.1029/2021WR030784

---

INPUT DATA

Here, we describe the basic input necessary to run D-CASCADE. Note that the Code is heavily customizable to adapt the model to each case study.
- ReachData: struct reporting for each reach of the network the attribute column variables. ReachData is structured in the same way as in the original CASCADE toolbox. Therefore, we suggest users read the documentation available for the original CASCADE toolbox, (https://github.com/mtangi/cascademodel), including the user manual, to learn more about the structure and field in ReachData and how it can be extracted from a digital elevation model (DEM). An example of ReachData is provided for both the Bega and 3S river case study.
- Q: matrix containing the discharge [m3/s] for each daily timestep in each reach of the network. An example of Q is provided for both the Bega and 3S river case study, for a 11-year simulation.
- psi: sediment classes defined in Krumbein phi [φ] scale (Krumbein and Sloss, 1963). Each sediment volume is classified as one of these classes. An example of how to define psi is shown in the scripts.
- Qbi_dep_in: cell structure reporting the initialized sediment volume in the deposit layer for the first timestep of the simulation, partitioned between the sediment classes. An example of Qbi_input is provided for the Bega river case study.
- Qbi_input: cell structure reporting the input sediment volume in each timestep for each reach from external sources (e.g., hillslope erosion or debris flows), partitioned between the sediment classes. An example of Qbi_input is provided for the 3S river case study.

---
OPTIONAL INPUT

D-CASCADE contains numerous customization options to allow the user to adapt the code to each case study. Different options may be specified as input to the main function using Name,Value pair arguments.

Transport capacity
As in the original CASCADE toolbox, the user can change the transport capacity formula used in the simulation to determine mobilized sediment volume in each reach, in each timestep. The available formulas are available as functions in folder transport_capacity_computation. The operating rule can be specified as a Name-Value pair argument to the primary function using the name 'tr_cap_equation' and as value, the ID of the function, listed below (see example in the script):
1) Parker and Klingeman (Parker and Klingeman, 1982)
2) Wilcock and Crowe (Wilcock and Crowe, 2003) 
3) Engelung and Hansen (Engelund and Hansen, 1967) (default)
4) Yang (Yang, 1984) 
5) Wong and Parker (Wong and Parker, 2006) 
6) Ackers and White (Ackers and White, 1973) 

For total sediment transport formulas (formulas 3 to 6), the user must specify a fractional computational method to partition the sediment transport between grain size classes. The chosen formula can be specified as a Name-Value pair argument to the primary function using the name 'partition_formula' and as value, the ID of the method, listed below (see example in the script):
1) Direct computation by the size fraction approach
2) BMF approach (Bed Material Fraction)
3) TCF approach (Transport Capacity Fraction) via the Molinas formula
These method are derived from Molinas and Wu (2000). 

The user can choose between two formulas to calculate the characteristic velocity of sediment. The chosen formula can be specified as a Name-Value pair argument to the primary function using the name 'velocity_formula' and as value, the ID of the method, listed below (see example in the script):

1) Virtual sediment velocity according to total transport capacity. This formulation estimates the sediment velocity via  the total transport capacity of the mobilized sediment volume, and guarantees an estimate of sediment velocity that is useful for large-scale modelling efforts (Hassan et al., 1991; Czuba and Foufoula-Georgiou, 2014; Czuba et al. ,2018). However, it provides a single estimation of sediment velocity for the mobilized sediment volume, without discriminating between different sediment classes.
2) Virtual sediment velocity by fractional transport capacity. This formula calculates sediment velocity by measuring transport capacity for each class independently from the grain size distribution of the mobilized material. Thus, interactions between grains of different classes do not influence sediment velocity. As a result, sediment velocity may vary between sediment classes.

---

DAM INPUT DATA

If reservoirs are present in the network, they can be integrated into the simulation. The model will dynamically simulate their release schedule, water and sediment storage, and the changes in the hydrology of the downstream channels. 
The dam position and reservoir features are defined in the DamData struct, reporting for each dam their main attributes:
- Name and Code contain the full and reduced name of the reservoir, respectively;
- portfolio reports "1" if the dam is to be included in the simulation, "0" otherwise;
- Node_ID and el_DEM contain the ID of the network node where the dam is located;
- design_discharge and min_discharge contain the minimum and design discharge of the dams. Releases higher than the design discharge will be partially delivered via the dam spillways;
- FSL contains the reservoir's full supply level at the start of the simulation. Sediment trapping may reduce this value;
- FSL_ResArea and FSL_ResVolume report the reservoir area and storage volume at FSL;
- WL_start contains the water level of the reservoir at the start of the simulation;
- WL_table is a matrix containing the reservoir area and storage volume at different supply levels. The model uses this table to derive reservoir level, area, and volume when one of these parameters is known. This is done via linear interpolation between the data reported in WL_table;
- mean_inflow reports the reservoir's mean water inflow [m3/s], calculated as an annual average. This parameter is used to deriver information on the reservoir trapping efficiency. 

The user can also specify what the management strategy of the reservoir is. This version of D-CASCADE supports three different management strategies that determine the reservoir release schedule. The operating rule can be specified as an Name-Value pair argument to the primary function using the name 'OrRule' and as value the ID of the operating rule, listed below (see example in the script):
1) Release equal to input: the dam releases exactly the input discharge, unless it is higher then the design discharge or lower than the minimum discharge;
2) Constant target: the model tries to keep the water level constant. Thus, if the reservoir level is in target, the release is equal to the estimated inflow, unless it exceeds the design discharge. The input needed is specified in field WL_constant in struct ORparameters.
3) Reservoir release curve based on biannual targets: the reservoir target level changes each timestep, according to a 4 parameter rule curve. Each day the rule determines the target water level, based on linear interpolation between the maximum and minimum reservoir levels (H1 and H2) and their respective dates (T1 and T2). The target level, in turn, determines the target release for the timestep. The needed input are the current date for the timestep (dateQ), and the dates and values of the target maximum and minimum SL (specified in field WL_twotargets in struct ORparameters).
For cases 2 and 3, in each timestep the model needs an estimate of the reservoir inflow to determine the release. This input is necessary since the decision-maker simulated in the model cannot know the daily inflow when it chooses the release, since the decision must be taken at the beginning of the timestep. In this version of D-CASCADE, the decision-maker assumes the inflow to the reservoir in each timestep to equal the inflow in the previous timestep.


---

OUTPUT DATA

D-CASCADE returns its outputs in the form of cell structure containing in the first column the name of the output and in the second the output itself.
In total, D-CASCADE returns two of these cell structures, plus one additional one if reservoir are added in the network:

DATA_OUTPUT contains the main aggregated output matrices obtained by postprocessing of the raw model outputs inside the main function. All data in data_output are txn matrices reporting the specified values for each timestep t and each reach n:
- Channel Width [m] contains the channel width. In this version of D-CASCADE, if no reservoir is present in the network, the channel width will not change over time, except in the flooded reaches upstream of the dams, where it will change according to the stored volume of the reservoir;
- Mobilized volume [m3] contains the total sediment volume mobilized, which depends on the reach transport capacity and the sediment availability in each timestep;
- Total sed in the reach [m3] contains the total sediment volume in the reach, which combines both the material in the deposit layer and the mobilized material;
- Reach Slope [-] contains the channel gradient. In this version of D-CASCADE, the channel gradient dynamically adapt to the changes in sediment storage in the deposit layer, using the methodology described in Czuba et al., 2017 (https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1002/2016JF003965);
- Deposit layer D50 [mm] contains the D50 of the deposit layer;
- Active layer D50 [mm] contains the D50 of the active layer;
- Daily transport capacity [m^3/day] reports the theoretical daily transport capacity, which determines the maximum amount of material that can be mobilized in each reach in a daily timestep;
- Deposited volume [m^3] contains the total volume of the material in the deposit layer;
- Discharge [m^3/s] contains the input discharge matrix Q;
- Total sed in the reach - per class [m^3] is a cx1 cell structure, which contains the total sediment volume in the reach of each sediment class c. This value combines the material in the deposit layer and the mobilized material.

EXTENDED_OUTPUT contains all the raw D-CASCADE outputs, as they are used in the model. Compared to the data in data_outputs, these outputs contain all the information stored in D-CASCADE, including sediment provenance and class composition for each material delivered to each reach, the deposit layer structure and others. These raw data are very large structures that must be processed before being used in data visualization. For more information, refer to the D-CASCADE Code.

DAM_OUTPUT contains all the D-CASCADE outputs referring to the conditions and operations of the reservoirs added in the network.
- DamDatabase_active contains all the information in DamDatabase only for the dams included in the dam portfolio considered in the simulation;
- ORparameters_active contains all the information in ORparameters only for the dams included in the dam portfolio considered in the simulation;
- Reservoir Volume reports the water volume stored in the reservoir in each timestep [m3];
- Release reports the water released by each reservoir in each timestep [m3];
- Volume at FSL reports the theoretical volume that can be stored in the reservoir at full supply level (FSL) in each timestep [m3]. The sediment volume stored in the reservoir causes this value to decrease from the original FSL volume, meaning part of the reservoir's capacity is lost due to sediment trapping;
- Reservoir sediment storage contains the total sediment volume [m3] stored in the reservoirs in each timestep. This volume is composed of the aggregated sediment deposits in each flooded reaches.
- Reach Energy Slope [-] contains the energy slope for each reach in each timestep. In this version of D-CASCADE, the Energy Slope differs from the channel gradient in the timestep only in flooded reaches.
- Flooded reaches ID at FSL contains, for each dam, the ID of the flooded or partially flooded reaches when the reservoir is at full supply level.

---

References:

Ackers, P. and White, W. R. (1973). Sediment transport: new approach and analysis. Journal of the Hydraulics Division, 99(11):2041–2060.

Engelund, F. and Hansen, E. (1967). A monograph on sediment transport in alluvial streams. Technical University of Denmark 0stervoldgade 10, Copenhagen K.

Krumbein, W. C. and Sloss, L. L. (1963). Stratigraphy and sedimentation. Technical report.

Czuba, J.A., Foufoula-Georgiou, E., 2014. A network-based framework for identifying potential synchronizations and amplifications of sediment delivery in river basins. Water Resources Research 50, 3826–3851.

Czuba, J.A., 2018. A lagrangian framework for exploring complexities of mixed-size sediment transport in gravel-bedded river networks. Geomorphology 321, 146–152.

Hassan, M.A., Church, M., Schick, A.P., 1991. Distance of movement of coarse particles in gravel bed streams. Water Resources Research 27, 503–511.

Parker, G. and Klingeman, P. C. (1982). On why gravel bed streams are paved. Water Resources Research, 18(5):1409–1423.

Tangi, M., Bizzi, S., Fryirs, K., & Castelletti, A. (2022). A dynamic, network scale sediment (dis)connectivity model to reconstruct historical sediment transfer and river reach sediment budgets. Water Resources Research, 58, e2021WR030784. https://doi.org/10.1029/2021WR030784

Wilcock, P. R. and Crowe, J. C. (2003). Surface-based transport model for mixed-size sediment. Journal of Hydraulic Engineering, 129(2):120–128.

Wong, M. and Parker, G. (2006). Reanalysis and correction of bed-load relation of meyer-peter and müller using their own database. Journal of Hydraulic Engineering, 132(11):1159–1168.

Yang, C. T. (1984). Unit stream power equation for gravel. Journal of Hydraulic Engineering, 110(12):1783–1797.
