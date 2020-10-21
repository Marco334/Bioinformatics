import json
# MOLECOLES FILE FROM https://pubchem.ncbi.nlm.nih.gov/compound/190#section=3D-Conformer
#JSON to Dataframe 
# Extract chemical JSON molecular data to 2 dataframe :
# df_ATM: Atoms componing the molecule data
# df_LEGAMI: bonds betewwn atoms data 


SOURCE_PTH   = str('................\\...........\\')
SOURCE_FL_NM = 'ADENINE_Conformer3D_CID_190.json'
FILE_V       = str(SOURCE_PTH + SOURCE_FL_NM )
# Opening JSON file
try:
  f = open(FILE_V, 'r')
  data = json.load(f)
  print(type(data))
  #return list from json
except:
   print("ERROR LOADING FILE")

# LEVEL 2 -------
ids     = list()
atomss  = list()
bondss  = list()
coordss = list()
propss  = list()
counts  = list()
#----------------

try:
for element in data['PC_Compounds']:
    ids     = element['id'    ]
    atomss  = element['atoms' ]
    bondss  = element['bonds' ]
    coordss = element['coords']
    propss  = element['props' ]
    counts  = element['count' ]
except:
    print("ERROR exploring JSON LEVEL 2 ")
   
try:
    df_ATM = DataFrame ([atomss["aid"],atomss["element"],coordss[0]['conformers'    ][0]['x'],coordss[0]['conformers'    ][0]['y'],coordss[0]['conformers'    ][0]['z'] ]).transpose()
    df_ATM.columns = ['ATM_INDX','ATM_TYPE','CD_X','CD_Y','CD_Z']
    df_ATM
except:
    print("ERROR Creating ATOM dataframe")
   
try:
    df_LEGAMI = DataFrame ([ bondss["aid1"],bondss["aid2"]   ]).transpose()
    df_LEGAMI.columns = ['BONDS_1','BONDS_2']
    df_LEGAMI
except:
    print("ERROR Creating BONDS dataframe")
