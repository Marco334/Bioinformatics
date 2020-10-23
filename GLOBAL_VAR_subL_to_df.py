import pandas

#VIEW and CONVERT subporsion of global variables starting with specific string in a dataframe
def slicedict(d, s):
    ''' filter variables name'''
    return {k:v for k,v in tt.items() if str(k).startswith(s) }
    
def FILT_GLOB_VARDICT_to_DF(tt, s):
   ''' filter global variables dictionary and create dataFrame'''
    data_GV=slicedict(tt, s)
    GV_df = pd.DataFrame.from_dict(data_GV, orient='index' )
    return(GV_df) 
    
tt = vars()
s = 'SUBSTRING_'
print(FILT_GLOB_VARDICT_to_DF(tt, s))
