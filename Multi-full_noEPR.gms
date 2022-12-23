*-------------------------------------------------------------------------------
*                                Defining the sets
*-------------------------------------------------------------------------------
sets
         h                                               /0*8759/
         first(h)        'first hour'
         last(h)         'last hour'
         m               'month'                         /1*12/
         w               'week'                          /1*52/
         d               'day'                           /1*365/
         tec             'all the technologies'          /offshore_f, offshore_g, onshore, PV, river, lake,nuclear, PHS, battery, OCGT, CCGT, CCGT-CCS, methanization, pyrogaseification, ngas, gastank, methanation, electrolysis, ITES, CTES, resistive, hpc, hpd, boilerc, boilerd, heat-net, EV_heavy, EV_light, EV_bus, EV_train, ICE_heavy, ICE_light, ICE_bus, hydrogen, H2_OCGT/
         gen(tec)        'the generation technologies'   /offshore_f, offshore_g, onshore, PV, river, lake,nuclear, methanization, pyrogaseification, ngas/
         elec(tec)       'electricity technologies'      /offshore_f, offshore_g, onshore, PV, river, lake,nuclear, PHS, battery, OCGT, CCGT, CCGT-CCS, H2_OCGT/
         gas(tec)        'gas technologiess'             /methanization, pyrogaseification, ngas, methanation, electrolysis, gastank/
         heat(tec)       'heat technologies'             /ITES, CTES, resistive, hpc, hpd, boilerc, boilerd/
         lowheat(heat)   'low temperature heat'          /ITES, CTES, resistive, hpc, hpd, boilerc, boilerd/
         midheat(heat)   'mid temperature heat'          /resistive, boilerd/
         highheat(heat)  'high temperature heat'         /boilerd/
         transport(tec)  'transportation technologies'   /EV_heavy,EV_light,EV_bus,EV_train,ICE_heavy,ICE_light,ICE_bus/
         ev(transport)   'the electric transportation'   /EV_heavy,EV_light,EV_bus,EV_train/
         ice(transport)  'the ICE transportation'        /ICE_heavy,ICE_light,ICE_bus/
         transporttype   'transport demand types'        /heavy, light, bus, train/
         elec_gen(tec)   'electricity generation tecs'   /offshore_f, offshore_g, onshore, PV, river, lake, nuclear/
         gas_gen(tec)    'gas generation technologies'   /methanization, pyrogaseification, ngas/
         str(tec)        'storage technologies'          /battery, PHS, gastank, hydrogen, ITES, CTES/
         str_elec(str)   'electric storage technologies' /battery, PHS/
         str_gas(str)    'gas storage technologies'      /gastank/
         str_heat(str)   'heat storage technologies'     /ITES, CTES/
         conv(tec)       'vector change technologies'    /OCGT, CCGT, CCGT-CCS, methanation, electrolysis, resistive, hpc, hpd, boilerc, boilerd, H2_OCGT/
         HP(conv)        'heat pumps'                    /hpc, hpd/
         from_elec(conv) 'conversion from electricity'   /methanation, electrolysis, hpc, hpd, resistive/
         from_gas(conv)  'conversion from gas'           /OCGT, CCGT, CCGT-CCS, boilerc, boilerd/
         vre(tec)        'variable tecs'                 /offshore_f, offshore_g, onshore, PV, river/
         frr(tec)        'technologies for upward FRR'   /lake, PHS, battery, OCGT, CCGT, CCGT-CCS, nuclear, H2_OCGT/
         central(tec)    'central heating technologies'  /hpc, boilerc, CTES/
         scenCO2         'CO2 tax scenarios'             /0, 100, 200, 300, 400, 500/
         scenVRE         'VRE cost scenarios'            /cheap, central, expensive/
         scenEPR         'nuclear power cost scenarios'  /cheap, central, expensive/
         scenBIO         'biogas options cost scenarios' /cheap, central, expensive/
         scenP2G         'methanation cost scenarios'    /cheap, central, expensive/
         costtype        'investment / operational costs'/capex,opex/
         biogas(tec)     'renewable gas options'         /methanization, pyrogaseification/
         type            'heat demand type'              /lowT, midT, highT/
;
first(h) = ord(h)=1;
last(h) = ord(h)=card(h);
alias(h,hh);
*-------------------------------------------------------------------------------
*                                Inputs
*-------------------------------------------------------------------------------
$offlisting
$ondelim
parameter day(h)
/
$include inputs/day.csv
/;
parameter week(h)
*last week is 9 days
/
$include inputs/week.csv
/;
parameter month(h) 'Production profiles of VRE'
/
$include  inputs/month.csv
/;
parameter load_factor(vre,h) 'Production profiles of VRE'
/
$include  inputs/vre_profiles2006f.csv
/;
parameter profile(transport,h) 'the charging profile of each type of transport'
/
$include inputs/t_profiles.csv
/;
parameter demand_elec(h) 'electricity demand profile in each hour in GWh-e'
/
$include inputs/elec_demand_2050.csv
/;
Parameter lake_inflows(m) 'monthly lake inflows in GWh'
/
$include  inputs/lake2006.csv
/ ;
parameter epsilon(vre) 'additional FRR requirement for variable renewable energies because of forecast errors'
/
$include  inputs/reserve_requirementsf.csv
/ ;
parameter demand_heat(h,type) 'heat demand profile in each hour in GWh-th'
*resource: GRTgas + TIGF
/
$include inputs/heat_demand_2050.csv
/;
parameter demand_transport(transporttype,h)    'hourly transport demand in GWh-useful'
/
$include inputs/demand_transport2050.csv
/;
parameter capa_ex(tec) 'existing capacities of the technologies by December 2017 in GW'
*Resource: RTE
/
$include  inputs/existing_capas.csv
/ ;
parameter capa_max(tec) 'maximum capacities of the technologies by December 2017 in GW'
*Resource: ADEME
/
$include  inputs/max_capas.csv
/ ;
parameter capex(tec) 'annualized capex cost in M€/GW/year'
/
$include  inputs/annuities_full.csv
/ ;
parameter capex_en(str)  'annualized capex cost in M€/GWh/year'
/
$include inputs/str_annuities_full.csv
/ ;
parameter fOM(tec) 'annualized fixed operation and maintenance costs M€/GW'
/
$include  inputs/fO&M_full.csv
/ ;
Parameter vOM(tec) 'Variable operation and maintenance costs in M€/GWh'
/
$include  inputs/vO&M_full.csv
/ ;
parameter cop_hp(h) 'hourly coefficient of performance for heat pumps'
/
$include inputs/HP_COP.csv
/ ;
*Parameter costsVRE(vre,costtype,scenVRE) 'VRE costs depending on scenarios'
*/
*$include  inputs/costsVRE.csv
*/ ;
*Parameter costsBIO(biogas,scenBIO) 'renewable gas costs depending on the scenarios'
*/
*$include  inputs/costsBIO.csv
*/ ;
*Parameter costsEPR(costtype,scenEPR) 'nuclear costs depending on the cost scenario'
*/
*$include  inputs/costsEPR.csv
*/ ;
*Parameter costsP2G(costtype,scenP2G) 'nuclear costs depending on the cost scenario'
*/
*$include  inputs/costsP2G.csv
*/ ;
$offdelim
$Onlisting
parameter eta_in(str) 'charging efifciency of storage technologies' /PHS 0.9, battery 0.9, gastank 1, ITES 1, CTES 1, hydrogen 1/;
parameter eta_out(str) 'discharging efficiency of storage technolgoies' /PHS 0.9, battery 0.95, gastank 0.99, ITES 0.9, CTES 0.75, hydrogen 0.98/;
parameter eta_conv(conv) 'the vector conversion efficiency of vector change options' /OCGT 0.45, CCGT 0.63, CCGT-CCS 0.53, methanation 0.6352, electrolysis 0.8, hpc 3.5, hpd 3.5, resistive 0.9, boilerc 0.9, boilerd 0.9, H2_OCGT 0.4/;
parameter eta(transport) 'the conversion ratio of electricity/gas to kms' /EV_heavy 0.847458,EV_light 8,EV_bus 0.847458, EV_train 1,ICE_heavy 0.347584,ICE_light 1.97,ICE_bus 0.347584/;
scalar pump_capa 'pumping capacity in GW' /9.3/;
scalar max_phs 'maximum volume of energy can be stored in PHS reservoir in TWh' /0.18/;
scalar load_uncertainty 'uncertainty coefficient for hourly demand' /0.01/;
scalar delta 'load variation factor'     /0.1/;
scalar tank_max 'maximum gas storage volume in TWh'    /134.6/   ;
scalar max_metha 'maximum yearly energy can be produced from methanisation in TWh'       /152/;
scalar max_pyro 'maximum yearly energy can be produced from pyrogasification in TWh'       /77/ ;
scalar max_central 'maximal yearly heat that can be satisfied by central heating in TWh' /151.19/;
scalar cf_nuc 'maximum capacity factor of nuclear power plants' /0.90/;
scalar ramp_rate 'maximum ramp up/down rate for nuclear power plant' /0.50/;
*scalar cf_ccgt 'maximum capaity factor of CCGT plant for a year' /0.90/;
*scalar cf_ccs /0.85/;
parameter CO2_tax(scenCO2) 'CO2 tax for each CO2 tax scenario' /0 0,100 100,200 200,300 300,400 400,500 500/;
parameter ngas_plus(scenCO2) 'the increase in price of natural gas because of CO2 tax' /0 0,100 22.95,200 45.9,300 68.85,400 91.8,500 114.75/;
parameter ccs_minus(scenCO2) 'the subtract in the carbon tax coming from CCS' /0 0,100 -29.2, 200 -58.4, 300 -87.6, 400 -116.8, 500 -146/;
scalar EV_batt 'the price energy related capex of EV batteries' /12.64/;
scalar CO2potential 'the maximal CO2 that can be stored yearly'        /93/;
parameter vOM0(*) /ngas 0.02354, CCGT-CCS 0.00576/;
parameter demand_hydrogen(h) 'hydrogen demand profile for industry - flat profile';
demand_hydrogen(h) = 4.56621;
*parameter demand_hydrogen(h) 'hydrogen demand profile for industry (40TWh) + aerial (57TWh) and marine (17TWh) transport - flat profile' /13.0137/;
parameter capacity_ex(str);
capacity_ex('hydrogen') = 3000;
capacity_ex('phs') = 80.16
*-------------------------------------------------------------------------------
*                                Model
*-------------------------------------------------------------------------------
variables        GENE(tec,h)     'hourly energy generation in TWh'
                 CAPA(tec)       'overal yearly installed capacity in GW'
                 STORAGE(str,h)  'hourly electricity input of battery storage GW'
                 CHARGE(transport,h)
                 SOC(str,h)      'state of charge of storage option in TWh'
                 S(str)          'storage charging capacity in GW'
                 CONVERT(conv,h) 'vector conversion from electricity to gas and the opposite in GW'
                 CAPACITY(str)   'energy volume of storage technologies in GWh'
                 RSV(frr,h)      'required upward frequency restoration reserve in GW'
                 EV_VOL          'the volume of electric vehicles in GWh'
                 COST            'final investment cost in b€'
positive variables GENE(tec,h),CAPA(tec),STORAGE(str,h),SOC(str,h), S(str), CAPACITY(str),RSV(frr,h),CONVERT(conv,h) ;
equations        gene_vre                'variables renewable profiles generation'
                 gene_capa               'capacity and genration relation for technologies'
                 capa_frr                'capacity needed for the secondary reserve requirements'
                 storing_const           'the definition of stored energy in the storage options'
                 conversion              'the mechanism of conversion'
                 conv_hp
                 storage_const           'storage in the first hour is equal to the storage in the last hour'
                 lake_res                'constraint on water for lake reservoirs'
                 stored_cap              'maximum energy that is stored in storage units'
                 storage_capa            'hourly storage should not be more than discharging capacity'
*                 nuc_cf                  'the yearly capacity factor of nuclear power plants should not pass 80%'
*                 nuc_up                  'Nuclear power plant upward flexibility flexibility'
*                 nuc_down                'Nuclear power plant downward flexibility flexibility'
*                 ccgt_cf                 'the yearly capacity factor of CCGT'
*                 ccs_cf
                 methanisation_max       'maximum yearly possible energy output from methanisation'
                 pyrogasification_max    'maximum yearly possible energy output from pyrogasification'
                 methanation_max         'maximum green CO2 that can be foud to produce methanation'
                 reserves                'FRR requirement'
                 heat_network            'the central technologies need heat network while decentrals dont need it'
                 heat_network_c          'the central technologies need heat network while decentrals dont need it'
                 max_heat                'the maximal heat that can be satisfied by central heating'
                 charging_profiles       'the charging profiles of EVs and ICEs'
                 trans_charge            'the weekly demand of transport should be satisfied during the week'
*                 max_CO2                 'an assumtion to limit the maximal storable CO2 in MtCO2/year'
                 hydrogen_electrolysis
                 ev_volume               'the weekly energy stored in electric vehicle should not be more than the available volume'
                 adequacy_elec           'supply/demand relation for the electricity demand'
                 adequacy_gas            'supply/demand relation for the gas demand'
                 adequacy_heat_low       'supply/demand relation for the heat demand'
                 adequacy_heat_mid       'supply/demand relation for the heat demand'
                 adequacy_t_heavy        'supply/demand equilibrium for transport demand'
                 adequacy_t_light        'supply/demand equilibrium for transport demand'
                 adequacy_t_bus          'supply/demand equilibrium for transport demand'
                 obj                     'the final objective function which is COST';
gene_vre(vre,h)..                GENE(vre,h)                     =e=     CAPA(vre)*load_factor(vre,h);
gene_capa(tec,h)..               CAPA(tec)                       =g=     GENE(tec,h);
capa_frr(frr,h)..                CAPA(frr)                       =g=     GENE(frr,h) + RSV(frr,h);
storing_const(h,h+1,str)..       SOC(str,h+1)                    =e=     SOC(str,h) + STORAGE(str,h)*eta_in(str) - GENE(str,h)/eta_out(str);
storage_const(str,first,last)..  SOC(str,first)                  =e=     SOC(str,last)+ STORAGE(str,last)*eta_in(str) - GENE(str,last)/eta_out(str); ;
lake_res(m)..                    lake_inflows(m)                 =g=     sum(h$(month(h) = ord(m)),GENE('lake',h))/1000;
stored_cap(str,h)..              SOC(str,h)                      =l=     CAPACITY(str);
storage_capa(str,h)..            STORAGE(str,h)                  =l=     CAPA(str);
conversion(h,conv)$(not HP(conv))..  GENE(conv,h)                    =e=     CONVERT(conv,h)*eta_conv(conv);
conv_hp(h,HP)..                  GENE(HP,h)                      =e=     CONVERT(HP,h)*cop_hp(h);
*nuc_cf..                         sum(h,GENE('nuclear',h))        =l=     CAPA('nuclear')*cf_nuc*8760;
*nuc_up(h,h+1)..                  GENE('nuclear',h+1) + RSV('nuclear',h+1) =l= GENE('nuclear',h) + ramp_rate*CAPA('nuclear')   ;
*nuc_down(h,h+1)..                GENE('nuclear',h+1) =g= GENE('nuclear',h)*(1 - ramp_rate)   ;
*ccgt_cf..                        sum(h,GENE('CCGT',h))           =l=     CAPA('CCGT')*cf_ccgt*8760;
*ccs_cf..                         sum(h,GENE('CCGT-CCS',h))       =l=     CAPA('CCGT-CCS')*cf_ccs*8760;
methanisation_max..              sum(h,GENE('methanization',h))  =l=     max_metha*1000;
pyrogasification_max..           sum(h,GENE('pyrogaseification',h))=l=   max_pyro*1000;
methanation_max..                sum(h,CONVERT('methanation',h)) =l=     sum(h,GENE('methanization',h))*3/7;
reserves(h)..                    sum(frr, RSV(frr,h))            =g=     sum(vre,epsilon(vre)*CAPA(vre))+ demand_elec(h)*load_uncertainty*(1+delta);
heat_network(central,h)..        CAPA('heat-net')                =g=     CAPA(central);
heat_network_c(h)..              GENE('heat-net',h)              =e=     sum(central,GENE(central,h));
max_heat..                       sum(h,GENE('heat-net',h))       =l=     max_central*1000;
charging_profiles(transport,h).. GENE(transport,h)               =l=     CAPA(transport)*profile(transport,h);
trans_charge(transport,w)..      sum(h$(ord(w)=week(h)),GENE(transport,h)) =g= sum(h$(ord(w) = week(h)), CHARGE(transport,h));
*max_CO2..                        CO2potential                    =g=     sum(h,GENE('ccgt-ccs',h))*292.4/1000000 ;
hydrogen_electrolysis(h)..       GENE('electrolysis',h)          =e=     STORAGE('hydrogen',h)+ demand_hydrogen(h)- GENE('hydrogen',h)+CONVERT('H2_OCGT',h);
ev_volume(w)..                   sum((h,ev)$(ord(w)=week(h)),GENE(ev,h)) =l= EV_VOL;
adequacy_elec(h)..               sum(elec,GENE(elec,h))          =g=     demand_elec(h) + sum(str_elec,STORAGE(str_elec,h)) + sum(from_elec,CONVERT(from_elec,h))+sum(ev,GENE(ev,h));
adequacy_gas(h)..                sum(gas,GENE(gas,h))            =e=     sum(str_gas,STORAGE(str_gas,h)) + sum(from_gas,CONVERT(from_gas,h))+sum(ice,GENE(ice,h));
adequacy_heat_low(h)..           sum(lowheat,GENE(lowheat,h))    =g=     sum(type, demand_heat(h,type)) + sum(str_heat,STORAGE(str_heat,h));
adequacy_heat_mid(h)..           sum(midheat,GENE(midheat,h))    =g=     demand_heat(h,'midT')+demand_heat(h,'highT');
adequacy_t_heavy(h)..            CHARGE('EV_heavy',h)*eta('EV_heavy')+CHARGE('ICE_heavy',h)*eta('ICE_heavy') =e= demand_transport('heavy',h);
adequacy_t_light(h)..            CHARGE('EV_light',h)*eta('EV_light')+CHARGE('ICE_light',h)*eta('ICE_light') =e= demand_transport('light',h);
adequacy_t_bus(h)..              CHARGE('EV_bus',h)*eta('EV_bus')+CHARGE('ICE_bus',h)*eta('ICE_bus') =e= demand_transport('bus',h);
obj..                            COST                            =e=     (sum(tec,(CAPA(tec)-capa_ex(tec))*capex(tec))+ sum(str,(CAPACITY(str)-capacity_ex(str))*capex_en(str))+sum(tec,(CAPA(tec)*fOM(tec))) +sum((tec,h),GENE(tec,h)*vOM(tec))+EV_VOL*EV_batt)/1000;
*-------------------------------------------------------------------------------
*                                Initial and fixed values
*-------------------------------------------------------------------------------
CAPA.fx('phs') = pump_capa;
CAPA.fx('river')= capa_ex('river');
CAPA.fx('lake') = 12.855;
S.fx('phs') = pump_capa;
CAPACITY.fx('phs') = max_phs*1000;
CAPA.up('offshore_f') = capa_max('offshore_f');
CAPA.up('offshore_g') = capa_max('offshore_g');
CAPA.up('onshore') = capa_max('onshore');
CAPA.up('pv') = capa_max('pv');
CAPACITY.fx('gastank') = tank_max*1000;
*two main equations defined here
CHARGE.fx('EV_train',h) = demand_transport('train',h)/eta('EV_train');
GENE.lo('boilerd',h) = demand_heat(h,'highT');
CAPACITY.lo('hydrogen') = capacity_ex('hydrogen');
*-----------------
CAPA.fx('nuclear')=0;
GENE.fx('nuclear',h)=0;
*-------------------------------------------------------------------------------
*                                Model options
*-------------------------------------------------------------------------------
model EOLES /all/;
option solvelink=0;
option RESLIM = 1000000;
option lp=CPLEX;
option Savepoint=1;
option solveopt = replace;
option limcol = 0;
option limrow = 0;
option SOLPRINT = OFF;
option solvelink=0;
$onecho > cplex.opt
objrng all
rhsrng all
rngrestart rng.inc
$offecho
EOLES.optfile=1; EOLES.dictfile=2;
*-------------------------------------------------------------------------------
*                                Solve statement
*-------------------------------------------------------------------------------
*$If exist EOLES_p.gdx execute_loadpoint 'EOLES_p';
parameter sumgene_elec; parameter sumgene_gas; parameter gene_tec(tec);
parameter dem_h2; parameter dem_elec; parameter dem_heat_high;parameter dem_heat_mid; parameter dem_heat_low; parameter dem_transport;
parameter elec_gas; parameter gas_elec; parameter elec_heat; parameter gas_heat; parameter elec_h2; parameter h2_elec;
parameter elec_transport; parameter gas_transport;
parameter str_loss_elec 'yearly power storage related loss in % of power production';
parameter conv_loss_elec 'overall energy loss from P2G conversion in % of power production';
parameter conv_loss_gas 'overall energy loss from G2P conversion in % of gas production';
*parameter lcoe(gen) 'levelized cost of energy production for each production technology in €/MWh';
*parameter lcos(str) 'levelized cost of energy storage for each storage technology in €/MWh';
*parameter lcoc(conv) 'levelized cost of conversion for each of vector change technologies in €/MWh';
parameter elec_spot_price(h) 'marginal cost'; parameter gas_spot_price(h) 'marginal cost';
parameter heat_price_high(h) 'marginal cost';parameter heat_price_mid(h) 'marginal cost';parameter heat_price_low(h) 'marginal cost';
parameter h2_spot_price(h) 'marginal cost';
parameter t_heavy_spot_price(h); parameter t_light_spot_price(h);parameter t_bus_spot_price(h);
parameter elec_marginal_cost 'average value over the year of power spot price in €/MWh-e';
parameter gas_marginal_cost 'average value over the year of gas spot price in €/MWh-th';
parameter heat_marginal_high 'average value over the year of gas spot price in €/MWh-th';
parameter heat_marginal_mid 'average value over the year of gas spot price in €/MWh-th';
parameter heat_marginal_low 'average value over the year of gas spot price in €/MWh-th';
parameter transport_marginal_cost; parameter h2_marginal_cost;
parameter G2P_bought 'overall yearly money spent to the gas bought for injection in the electricity in M€';
parameter P2G_bought 'overall yearly money spent to the electricity bought for injection in the gas in M€';
parameter G2H_bought; parameter P2H_bought; parameter P2T_bought; parameter G2T_bought;
parameter lcoe_elec; parameter lcoe_gas; parameter lcoe_heat; parameter lcoe_h2;
parameter lc_elec;
*parameter all_losses 'load curtailment, storage losses and conversion losses';
*parameter cf(tec) 'load factor of all technologies';
parameter net_emission 'in MtCO2/year'; parameter captured_CO2 'in MtCO2/year'; parameter stored_CO2;
parameter technical_cost;
file hourly_profiles /'outputs/hourly_profiles_noNUC_n.csv'/ ;
file results_csv /'outputs/summary_noNUC_n.csv'/;
parameter nSTORAGE(str,h);
parameter nCONVERT(conv,h);
put hourly_profiles;
put 'CO2_scen','hour'; loop(tec, put tec.tl;) put 'demand_elec', 'ElecStr','Pump','H2_str','demand_hydrogen','gastank','demand_heat','CTES','ITES','elec_market','gas_market','heat_market_high','heat_market_mid','heat_market_low',loop(conv, put conv.tl;) put 'OK'/ ;
put results_csv;
results_csv.pc=5;
results_csv.pw=32767;
*put 'scenCO2','scenVRE','scenEPR','scenBIO','scenP2G','cost'; loop(tec, put tec.tl;) loop(tec,put tec.tl;) put 'LCOE-e','LCOE-g','spot-e','spot-gas','LC_elec','str_loss_elec','conv_loss_P2G','conv_loss_G2P','net_emission','CO2_captured'/;
put 'scenCO2','cost', 'technical cost'; loop(tec, put tec.tl;) loop(tec,put tec.tl;)put 'EV_VOL'; put 'LCOE-e','LCOE-g','LCOE-h','spot-e','spot-gas','spot-heat','net_emission','CO2_captured'/;
loop(scenCO2,
vOM('ngas') = vOM0('ngas')+ngas_plus(scenCO2)/1000;
vOM('CCGT-CCS') = vOM0('CCGT-CCS')+ccs_minus(scenCO2)/1000;
*capex('offshore') = costsVRE('offshore','capex',scenVRE);
*fOM('offshore') = costsVRE('offshore','opex',scenVRE);
*capex('onshore') = costsVRE('onshore','capex',scenVRE);
*fOM('onshore') = costsVRE('onshore','opex',scenVRE);
*capex('pv') = costsVRE('pv','capex',scenVRE);
*fOM('pv') = costsVRE('pv','opex',scenVRE);
*capex('nuclear') = costsEPR('capex',scenEPR);
*fOM('nuclear') = costsEPR('opex',scenEPR);
*vOM('methanisation') = costsBIO('methanisation',scenBIO);
*vOM('pyrogasification') = costsBIO('pyrogasification',scenBIO);
*capex('methanation') = costsP2G('capex',scenP2G);
*fOM('methanation') = costsP2G('opex',scenP2G);
Solve EOLES using lp minimizing COST;
sumgene_elec = sum((elec_gen,h),GENE.l(elec_gen,h))/1000;
sumgene_gas = sum((gas_gen,h),GENE.l(gas_gen,h))/1000;
gene_tec(tec) = sum(h,GENE.l(tec,h))/1000;
dem_h2 = sum(h,demand_hydrogen(h))/1000;
dem_elec = sum(h,demand_elec(h))/1000;
dem_heat_high = sum(h,demand_heat(h,'highT'))/1000;
dem_heat_mid = sum(h,demand_heat(h,'midT'))/1000;
dem_heat_low = sum(h,demand_heat(h,'lowT'))/1000;
dem_transport = sum((transporttype,h),demand_transport(transporttype,h))/1000;
elec_gas = sum(h,CONVERT.l('methanation',h))/1000;
elec_h2 = sum(h,CONVERT.l('electrolysis',h))/1000;
gas_elec = sum(h,(CONVERT.l('OCGT',h)+CONVERT.l('CCGT-CCS',h)+CONVERT.l('CCGT',h)))/1000;
h2_elec = sum(h,CONVERT.l('H2_OCGT',h))/1000;
elec_heat = sum(h,(CONVERT.l('resistive',h)+CONVERT.l('hpc',h)+CONVERT.l('hpd',h)))/1000;
gas_heat = sum(h,(CONVERT.l('boilerd',h)+CONVERT.l('boilerc',h)))/1000;
elec_transport = sum((h,ev),CHARGE.l(ev,h))/1000;
gas_transport = sum((h,ice),CHARGE.l(ice,h))/1000;
str_loss_elec = (sum((str_elec,h),STORAGE.l(str_elec,h))-sum(str_elec,gene_tec(str_elec)*1000))*100/((sumgene_elec-elec_gas - elec_heat - elec_transport )*1000);
*lcoe(gen) = (CAPA.l(gen)*(fOM(gen)+capex(gen))+ (gene_tec(gen)*vOM(gen)*1000))/gene_tec(gen);
*lcos(str) = (CAPA.l(str)*(fOM(str)+capex(str))+ (gene_tec(str)*vOM(str)*1000) + CAPACITY.l(str)*capex_en(str))/gene_tec(str);
*lcoc(conv) = (CAPA.l(conv)*(fOM(conv)+capex(conv))+ (gene_tec(conv)*vOM(conv)*1000))/gene_tec(conv);
elec_spot_price(h) = 1000000*adequacy_elec.m(h);
gas_spot_price(h) = 1000000*adequacy_gas.m(h);
heat_price_high(h) = 1000000*GENE.m('boilerd',h);
heat_price_mid(h) = 1000000*adequacy_heat_mid.m(h);
heat_price_low(h) = 1000000*adequacy_heat_low.m(h);
t_heavy_spot_price(h) = 1000000*adequacy_t_heavy.m(h);
t_light_spot_price(h) = 1000000*adequacy_t_light.m(h);
t_bus_spot_price(h) = 1000000*adequacy_t_bus.m(h);
h2_spot_price(h) = 1000000*hydrogen_electrolysis.m(h);
elec_marginal_cost = sum(h,elec_spot_price(h))/8760;
gas_marginal_cost = sum(h,gas_spot_price(h))/8760;
heat_marginal_high = sum(h,heat_price_high(h))/8760;
heat_marginal_mid = sum(h,heat_price_mid(h))/8760;
heat_marginal_low = sum(h,heat_price_low(h))/8760;
h2_marginal_cost = sum(h,h2_spot_price(h))/8760;
G2P_bought = sum(h,gas_spot_price(h)*(CONVERT.l('OCGT',h)+CONVERT.l('CCGT-CCS',h)+CONVERT.l('CCGT',h)))/1000;
P2G_bought = sum(h,elec_spot_price(h)*(CONVERT.l('methanation',h)+CONVERT.l('electrolysis',h)))/1000;
G2H_bought = sum(h,gas_spot_price(h)*(CONVERT.l('boilerc',h))+CONVERT.l('boilerd',h))/1000;
P2H_bought = sum(h,elec_spot_price(h)*(CONVERT.l('hpc',h)+CONVERT.l('hpd',h)+CONVERT.l('resistive',h)))/1000;
P2T_bought = sum((h,ev),elec_spot_price(h)*GENE.l(EV,h))/1000;
G2T_bought = sum((h,ice),gas_spot_price(h)*GENE.l(ice,h))/1000;
lcoe_elec =( sum(elec,(CAPA.l(elec)*(fOM(elec)+capex(elec)))+ (gene_tec(elec)*vOM(elec)*1000)) + sum(str_elec, CAPACITY.l(str_elec)*capex_en(str_elec)) +G2P_bought)/(sumgene_elec+gene_tec('OCGT')+gene_tec('CCGT')+gene_tec('CCGT-CCS')) ;
lcoe_gas = (sum(gas,(CAPA.l(gas)*(fOM(gas)+capex(gas)))+ (gene_tec(gas)*vOM(gas)*1000))+P2G_bought)/(sumgene_gas+gene_tec('methanation')+gene_tec('electrolysis'));
lcoe_heat = (sum(heat,(CAPA.l(heat)*(fOM(heat)+capex(heat)))+ (gene_tec(heat)*vOM(heat)*1000))+P2H_bought+G2H_bought)/(gene_tec('resistive')+gene_tec('hpc')+gene_tec('hpd')+gene_tec('boilerc')+gene_tec('boilerd'));
lc_elec = (sumgene_elec - elec_gas - elec_heat - elec_transport - dem_elec)*100/(sumgene_elec - elec_gas - elec_heat - elec_transport) - str_loss_elec;
*all_losses = (sumgene_elec + sumgene_gas - dem_elec - dem_gas)*100/(sumgene_elec+sumgene_gas);
*cf(tec) = gene_tec(tec)*1000/(8760*CAPA.l(tec));
net_emission = 0.2295*gene_tec('ngas')-0.292*gene_tec('CCGT-CCS');
captured_CO2 = 0.292*gene_tec('CCGT-CCS');
technical_cost = COST.l + (-gene_tec('ngas')*0.2295 + gene_tec('CCGT-CCS')*0.292)*CO2_tax(scenCO2)/1000;
nSTORAGE(str,h) = 0 - STORAGE.l(str,h);
nCONVERT(conv,h) = 0 - CONVERT.l(conv,h);
hourly_profiles.pc=5;
display cost.l; display capa.l; display sumgene_elec; display sumgene_gas;
display gene_tec; display dem_elec; display dem_h2; display dem_heat_high; display dem_heat_mid; display dem_heat_low; display dem_transport;
display capacity.l; display elec_gas; display gas_elec; display elec_heat; display gas_heat;
*display str_loss_elec; display conv_loss_elec; display conv_loss_gas;
display lcoe_elec;  display lcoe_gas; display lcoe_heat;
display net_emission; display captured_CO2;
*display lcoe; display lcos; display lcoc;
display lc_elec;
display elec_marginal_cost; display gas_marginal_cost; display heat_marginal_high, heat_marginal_mid, heat_marginal_low;
display G2P_bought; display P2G_bought; display G2H_bought; display P2H_bought; display P2T_bought; display G2T_bought;
*-------------------------------------------------------------------------------
*                                Outputs
*-------------------------------------------------------------------------------
put results_csv;
results_csv.pc=5;
*put scenCO2.tl,scenVRE.tl,scenEPR.tl,scenBIO.tl,scenP2G.tl, COST.l, loop(tec, put CAPA.l(tec);) loop(tec,put gene_tec(tec);) put lcoe_elec,lcoe_gas,elec_marginal_cost,gas_marginal_cost,lc_elec,str_loss_elec, conv_loss_elec,conv_loss_gas, net_emission, captured_CO2 /;
put scenCO2.tl, COST.l, technical_cost, loop(tec, put CAPA.l(tec);) loop(tec,put gene_tec(tec);)put EV_VOL.l; put lcoe_elec,lcoe_gas,lcoe_heat,elec_marginal_cost,gas_marginal_cost,heat_marginal_high,heat_marginal_mid,heat_marginal_low,net_emission, captured_CO2 /;
put hourly_profiles;
hourly_profiles.pc=5;
loop (h,
put scenCO2.tl, put h.tl; loop(tec, put GENE.l(tec,h);) put demand_elec(h); put nSTORAGE('battery',h), nSTORAGE('phs',h), nSTORAGE('hydrogen',h), demand_hydrogen(h); put nSTORAGE('gastank',h); put nSTORAGE('CTES',h),nSTORAGE('ITES',h),elec_spot_price(h),gas_spot_price(h), heat_price_high(h),heat_price_mid(h),heat_price_low(h); loop(conv, put nCONVERT(conv,h);) put 'OK'/
;);
);
*-------------------------------------------------------------------------------
*                                The End :D
*-------------------------------------------------------------------------------
