#!/bin/bash            
source activate basemap

python -u ./main_scale_to_radar_OFP.py LE_D1_1km_30sec_OFP_V2 > LE_D1_1km_30sec_OFP_V2.log  
python -u ./main_scale_to_radar_OFP.py LE_D1_1km_1min > LE_D1_1km_1min.log           
python -u ./main_scale_to_radar_OFP.py LE_D1_1km_1min_4D > LE_D1_1km_1min_4D.log        
python -u ./main_scale_to_radar_OFP.py LE_D1_1km_2min_OFP_V2 > LE_D1_1km_2min_OFP_V2.log    
python -u ./main_scale_to_radar_OFP.py LE_D1_1km_5min_OFP_V2 > LE_D1_1km_5min_OFP_V2.log    
python -u ./main_scale_to_radar_OFP.py LE_D1_1km_5min_4D_OFP_V2 > LE_D1_1km_5min_4D_OFP_V2.log 


