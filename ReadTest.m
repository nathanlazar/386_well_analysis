%% Read and adjust images
pref = 'C:\Users\lazar\Desktop\Thesis\CSMA-MEMA\386images\HCC1143.KRTVIM.NVS1.dmso\C5--W00053--P00001--Z00000--T00000--'
dapi = imadjust(imread(strcat(pref, 'DAPI.tif')));
alexa647 = imadjust(imread(strcat(pref, 'Alexa 647.tif')));
alexa488 = imadjust(imread(strcat(pref, 'Alexa 488.tif')));
alexa568 = imadjust(imread(strcat(pref, 'Alexa 568.tif')));