
import json
import pandas as pd
import mendeleev
import bpy
import os
import numpy as np
from periodictable import *
#MOLECOLE FROM https://pubchem.ncbi.nlm.nih.gov/compound/190#section=3D-Conformer

#ARTICOLO FIGATA ALLUCINANTE :https://medium.com/@behreajj/scripting-curves-in-blender-with-python-c487097efd13

#-------- info su tavola periodica
#https://stackoverflow.com/questions/20440343/how-to-provide-periodic-table-information-to-python-module
#https://docs.eyesopen.com/toolkits/python/oechemtk/periodictable.html

#------------------ CILINDRI
#https://sinestesia.co/blog/tutorials/python-tubes-cilinders/
#https://blender.stackexchange.com/questions/5898/how-can-i-create-a-cylinder-linking-two-points-with-python

#--------------- ARTICOLI SU PROGRAMMI SIMILI
#https://patrickfuller.github.io/molecules-from-smiles-molfiles-in-blender/
#--TO_DO:
# - Trovare atomi di contatto tra molecole
# - correzione rotazione su Y
# capire orientamento eliche
# - capire geometrie gradi rotazione eliche
# - Lettore sequena
# ESEMPIO DNA SEQUENCE 
#----https://www.bioinformatics.org/sms2/random_dna.html 


def LOAD_FILE(FILE_V):
                try:
                  f = open(FILE_V, 'r')
                  data = json.load(f)
                  print("MOLECULE FILE LOAD                           - DONE")
                  #return list from json
                except:
                   print("LOADING FILE - ERROR")
                return(data)

def CLLCT_GEN(NM):
                '''COLLECTION MANAGEMENT'''
                collection_T = bpy.data.collections.new("NUCLEOTIDE_"+ NM )
                bpy.context.scene.collection.children.link(collection_T)
                return(1)

def Glob_CLLCT_GEN():
                '''COLLECTION MANAGEMENT'''
                collection_T = bpy.data.collections.new("DS_Adenilato_0")
                bpy.context.scene.collection.children.link(collection_T)
                #-----------   --   --   --   --   --   --   --   --   --
                collection_T = bpy.data.collections.new("DS_Gualinato_0")
                bpy.context.scene.collection.children.link(collection_T)
                #-----------   --   --   --   --   --   --   --   --   --
                collection_T = bpy.data.collections.new("DS_Timidilato_0")
                bpy.context.scene.collection.children.link(collection_T)
                #-----------   --   --   --   --   --   --   --   --   --
                collection_T = bpy.data.collections.new("DS_Citidilato_0")
                bpy.context.scene.collection.children.link(collection_T)
                #-----------   --   --   --   --   --   --   --   --   --
                collection_T = bpy.data.collections.new("BASIC_PAIRS")
                bpy.context.scene.collection.children.link(collection_T)
                #-----------   --   --   --   --   --   --   --   --   --
                collection_T = bpy.data.collections.new("BASIC_PAIRS")
                bpy.context.scene.collection.children.link(collection_T)
                collection_T = bpy.data.collections.new('ATOMS')
                bpy.context.scene.collection.children.link(collection_T)
                collection_T = bpy.data.collections.new("BOUNDS")
                bpy.context.scene.collection.children.link(collection_T)
                collection_T = bpy.data.collections.new("DESCRIPTION_LAB")
                bpy.context.scene.collection.children.link(collection_T)
                return('COLLECTIONS READY')

def DF_GENERATOR(data):
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
                           df_ATM = pd.DataFrame([atomss["aid"],atomss["element"],coordss[0]['conformers'    ][0]['x'],coordss[0]['conformers'    ][0]['y'],coordss[0]['conformers'    ][0]['z'] ]).transpose()
                           df_ATM.columns = ['ATM_INDX','ATM_TYPE','CD_X','CD_Y','CD_Z']
                           #df_ATM
                       except:
                           print("ERROR Creating ATOM dataframe")
                       try:
                           df_LEGAMI = pd.DataFrame([ bondss["aid1"],bondss["aid2"]   ]).transpose()
                           df_LEGAMI.columns = ['BONDS_1','BONDS_2']
                           #df_LEGAMI
                       except:
                           print("ERROR Creating BONDS dataframe")
                       return(df_ATM,df_LEGAMI)


def MOLECULE_SFARE_GENERATOR(df_ATM,X_PLUS_FACTOR,COLLECT_NM):
                       for index, row in df_ATM.iterrows():
                           AT_INDEX   = int(row['ATM_INDX'])
                           AT_TYPE_NM = int(row['ATM_TYPE'])
                           global AT_RAD
                           global AT_NME
                           global AT_CLR
                           global AT_SYM
                           global AT_RGB
                           AT_RAD     = mendeleev.element( AT_TYPE_NM ).atomic_radius
                           AT_NME     = mendeleev.element( AT_TYPE_NM ).name
                           AT_CLR     = mendeleev.element( AT_TYPE_NM ).cpk_color
                           AT_SYM     = mendeleev.element( AT_TYPE_NM ).symbol
                           #print(str(AT_INDEX) + ' -------------------------------------------- '   )
                           #print(str(AT_RAD) + '  ' + str(AT_NME) + '  ' + str(AT_CLR) + '  ' + str(AT_SYM) )
                           AT_X       = row['CD_X']  #'CD_X']
                           AT_Y       = row['CD_Y']  #'CD_Y']
                           AT_Z       = row['CD_Z']  #'CD_Z']
                           global COLOR_NM
                           COLOR_NM   = "COLOR_" + str(int(AT_INDEX))
                           #-- COLORE ATOMO
                           #print(AT_CLR)
                           AT_RGB    = hex_to_rgb(AT_CLR)
                           script_C1 = COLOR_NM +" = bpy.data.materials.new(name='"+ COLOR_NM + "')"
                           exec(script_C1, globals())
                           #print(AT_RGB)
                           script_C2 = COLOR_NM +".diffuse_color = ( AT_RGB[0] , AT_RGB[1],AT_RGB[2] , 0 )"
                           exec(script_C2, globals())
                           #print('START BUILD SFARE')
                           ATOM_SFARE_GENERATOR(AT_X+X_PLUS_FACTOR,AT_Y,AT_Z,AT_RAD ,AT_NME ,AT_CLR ,AT_SYM,str(int(AT_INDEX)),COLOR_NM,COLLECT_NM)
                           del AT_RAD
                           del AT_NME
                           del AT_CLR
                           del AT_SYM
                           bpy.ops.object.select_all(action='DESELECT')
                       print('MOLECULE_SFARE_GENERATOR                     - DONE')
                       return()

def ATOM_SFARE_GENERATOR(AT_X,AT_Y,AT_Z,ATP_RAD,AT_NME ,AT_CLR ,AT_SYM,AT_INDEX,COLOR_NM,COLLECT_NM):
                       #print( ' RADIUS:' +str(ATP_RAD) + ' \n ' +  AT_NME + ' \n ' +str(AT_CLR) + ' \n ' +str(AT_SYM)+ ' \n ' +str(AT_INDEX)+ ' \n ' +str(COLOR_NM))
                       #COLOR_NM
                       ATP_RAD = ATP_RAD/100
                       bpy.ops.mesh.primitive_ico_sphere_add(radius= ATP_RAD, subdivisions=5 ,location=(AT_X,AT_Y,AT_Z))
                       obj = bpy.context.selected_objects[-1]
                       obj.name = str(AT_INDEX) +"_"+ str(COLLECT_NM)[3:5] + "_ATOMO_" + str(AT_NME)
                       #bpy.context.active_object.data.materials.append(COLOR_NM)
                       script_1 = str("bpy.context.active_object.data.materials.append("+COLOR_NM +")")
                       exec(script_1)
                       bpy.data.collections[COLLECT_NM].objects.link(obj)
                       bpy.ops.object.select_all(action='DESELECT')
                       return(str(AT_INDEX) + '  '+ str(AT_NME) +'  HAVE BEEN CREATED' )

def hex_to_rgb(color_str):
                       # supports '123456', '#123456' and '0x123456'
                       (r,g,b), a = map(lambda component: component / 255, bytes.fromhex(color_str[-6:])), 1.0
                       return (r,g,b,a)

def PREPARE_DF_BOUNDS(df_ATM,df_LEGAMI):
                       '''Double INNER JOIN BTW Dataframes to rretrive bound start end coordinates '''
                       ff = df_LEGAMI.merge(df_ATM , left_on='BONDS_1', right_on='ATM_INDX', how='inner')[['BONDS_1', 'CD_X','CD_Y','CD_Z'  , 'BONDS_2']]
                       ff.rename(columns = {'BONDS_1':'BONDS_1', 'CD_X':'A1_CD_X','CD_Y':'A1_CD_Y','CD_Z':'A1_CD_Z' , 'BONDS_2':'BONDS_2'}, inplace = True)
                       fg = ff.merge(df_ATM , left_on='BONDS_2', right_on='ATM_INDX', how='inner')[['BONDS_1','A1_CD_X', 'A1_CD_Y', 'A1_CD_Z' ,  'BONDS_2', 'CD_X','CD_Y','CD_Z' ]]
                       fg.rename(columns = {'BONDS_1':'BONDS_1', 'CD_X':'A1_CD_X','CD_Y':'A1_CD_Y','CD_Z':'A1_CD_Z' , 'BONDS_2':'BONDS_2', 'CD_X':'A2_CD_X','CD_Y':'A2_CD_Y','CD_Z':'A2_CD_Z'}, inplace = True)
                       print('DF BOUNDS PREP                               - DONE')
                       return(fg)

def CREATE_BOUNDS_1(df_BND,X_PLUS_FACTOR,COLLECT_NM):
                       global GRAY
                       for index, row in df_BND.iterrows():
                           bound_name = str(COLLECT_NM)[3:5]+'_BOUND_'+str(int(row['BONDS_1']))+'_' +str(int(row['BONDS_2']))
                           BNDS_EXTR  = [float(row['A1_CD_X'])+float(X_PLUS_FACTOR),float(row['A1_CD_Y']), float(row['A1_CD_Z']), float(row['A2_CD_X'])+float(X_PLUS_FACTOR),float(row['A2_CD_Y']), float(row['A2_CD_Z'])  ]
                           bpy.context.scene.cursor.location = (0,0,0)
                           bpy.ops.curve.primitive_bezier_curve_add()
                           obj = bpy.context.object
                           obj.data.dimensions       = '3D'
                           obj.data.fill_mode        = 'FULL'
                           obj.data.bevel_depth      = 0.1
                           obj.data.bevel_resolution = 0
                           # set first point to centre of sphere1
                           obj.data.splines[0].bezier_points[0].select_right_handle = True
                           obj.data.splines[0].bezier_points[0].co = (BNDS_EXTR[0],BNDS_EXTR[1],BNDS_EXTR[2])
                           obj.data.splines[0].bezier_points[0].handle_right_type = 'VECTOR'
                           # set second point to centre of sphere2
                           obj.data.splines[0].bezier_points[1].select_left_handle = True
                           obj.data.splines[0].bezier_points[1].co = (BNDS_EXTR[3],BNDS_EXTR[4],BNDS_EXTR[5])
                           obj.data.splines[0].bezier_points[1].handle_left_type = 'VECTOR'
                           obj = bpy.context.selected_objects[-1]
                           obj.name = bound_name
                           CURRENT  = bpy.data.objects.get(bound_name)
                           bpy.ops.object.convert(target='MESH', keep_original= False) #MASH CONVERSION
                           obj = bpy.context.active_object
                           bpy.context.active_object.data.materials.append( GRAY )
                           bpy.data.collections[COLLECT_NM].objects.link(obj)
                           bpy.ops.object.select_all(action='DESELECT')
                       print ('CREATE_BOUNDS                                - DONE')
                       bpy.ops.object.select_all(action='DESELECT')
                       return()

def CREATE_LABELS(X_PLUS_FACTOR,MOLECULE_NM_L,index_l):
                       #print('CREATE_LABELS - STARTED')
                       '''NATIONS LABEL MANAGEMENT'''
                       global TWD_CN
                       bpy.ops.object.select_all(action='DESELECT')
                       font_curve = bpy.data.curves.new(type="FONT",name= MOLECULE_NM_L  )
                       font_curve.body = str(str(MOLECULE_NM_L)).strip()    #contenuto
                       font_obj = bpy.data.objects.new( MOLECULE_NM_L, font_curve)
                       font_obj.location = (X_PLUS_FACTOR,-10, 0 )
                       font_obj.name = str(  "LABEL_" + MOLECULE_NM_L  )
                       bpy.data.collections['DESCRIPTION_LAB'].objects.link(font_obj)
                       #bpy.ops.object.select_all(action='DESELECT')
                       print('CREATE_LABELS                                - DONE')
                       return()

def BASIC_MOLECULE_JOIN(CLLCT_LIST):
                       for i in CLLCT_LIST:
                            COLL_NM = i
                            COLL_NM
                            for obj in bpy.data.collections[COLL_NM].all_objects:
                               if obj.type == 'MESH':
                                  obj.select_set( state = True, view_layer = None)
                                  bpy.context.view_layer.objects.active = obj
                                  #bpy.context.scene.objects.active = obj
                               else:
                                  obj.select_set( state = False, view_layer = None)
                            bpy.ops.object.join()
                            obj.name = "M0" + str(COLL_NM)[2:-2]
                            bpy.ops.object.select_all(action = 'DESELECT')
                       print('BASIC_MOLECULE_JOIN                          - DONE')
                       return()

def UNLINK_BASIC_MOLECULE_CL(CLLCT_LIST):
                       for i in CLLCT_LIST:
                            COLL_NM = i
                            for obj in bpy.data.collections[COLL_NM].all_objects:
                                  obj.select_set( state = True, view_layer = None)
                                  bpy.context.view_layer.objects.active = obj
                                  bpy.data.collections['Collection'].objects.unlink(obj)
                            bpy.ops.object.select_all(action = 'DESELECT')
                       print('Molecules M0 unlinked from Collection        - DONE')
                       return()

def UNLINK_BASIC_MOLECULE_GN(CLLCT_LIST):
                       for i in CLLCT_LIST:
                            COLL_NM = i
                            for obj in bpy.data.collections[COLL_NM].all_objects:
                                  obj.select_set( state = True, view_layer = None)
                                  bpy.context.view_layer.objects.active = obj
                                  bpy.context.scene.collection.objects.unlink(obj)
                            bpy.ops.object.select_all(action = 'DESELECT')
                       print('Molecules M0 unlinked from Master collection - DONE')
                       return()

def OBJ_MASS_C_MOVE_G():
                       scene       = bpy.context.scene
                       #global obj_M0_name
                       #for NM in obj_M0_name:
                       for obj in bpy.context.scene.objects:
                            #bpy.ops.object.select_pattern(pattern = NM )
                            bpy.ops.object.select_pattern(pattern = obj.name )
                            bpy.ops.object.origin_set(type   = 'ORIGIN_CENTER_OF_VOLUME')
                            bpy.ops.object.origin_set(type   = 'ORIGIN_CENTER_OF_MASS')
                            bpy.ops.object.select_all(action = 'DESELECT')
                       bpy.ops.object.select_all(action = 'DESELECT')
                       print('Molecules M0 MASS CENTER ROTATION            - DONE')
                       return()
 
def OBJ_ROTATE_G():
                       scene      = bpy.context.scene
                       obj_g_name = "M0_Gualinato"
                       template_ob = bpy.data.objects[obj_g_name].select_set(True)
                       obj = bpy.context.selected_objects[-1]
                       obj.rotation_euler.z = radians(180)
                       bpy.ops.object.transform_apply( rotation = True )
                       bpy.ops.object.select_all(action = 'DESELECT')
                       print('OBJ_ROTATE_G                                 - DONE')
                       return()
 
def OBJ_ROT_CORRRECTION():
                       ''' ROTATION CORRECTION BEFORE Pairs creation '''
                       bpy.ops.object.select_all(action = 'DESELECT')
                       #angle = atan2(y1 - y2, x1 - x2)
                       #angle
                       #angle * 180 / math.pi
                       for NM in obj_M0_name:
                           template_ob = bpy.data.objects[NM].select_set(True)
                           obj = bpy.context.selected_objects[-1]
                           obj.rotation_euler.x =  - 0.18541805585921572
                           #bpy.ops.object.transform_apply( rotation = True )
                           bpy.ops.object.select_all(action = 'DESELECT')
                       print('OBJ_ROT_COR                                  - DONE')
                       return()

def PROCESS_ROADMAP_PREP(FUNCTION_LST):
                       PROS_DF = pd.DataFrame(FUNCTION_LST, columns=['PROCESS'])
                       PROS_DF["EXCT_FL"] = np.nan
                       PROS_DF
                       print('PROCESS_ROADMAP_PREP                         - DONE')
                       return(PROS_DF)

def COLUMN_OBJ():
                       global obj_M0_name
                       for NM in obj_M0_name:
                           obj_M0_name.index(NM)
                           bpy.data.objects[NM].location = (0,(0+obj_M0_name.index(NM)*10),0)
                           bpy.data.objects[NM].rotation_euler.y = radians(180)
                       print('MOVE_IN_COLUMN_M0_OBJ                        - DONE')
                       return()

def CREATE_PAIRS():
                       global obj_M0_name
                       global PAIRS_DICT
                       AIRS_DICT = {}
                       bpy.ops.object.select_all(action = 'DESELECT')
                       for NM in obj_M0_name:
                           #NM_R = bpy.data.objects[NM].name
                           OR_OBJ_LC = bpy.data.objects[NM].location
                           if NM == 'M0_Adenilato':
                              COMPARE = 'M0_Timidilato'
                           elif NM == 'M0_Gualinato':
                              COMPARE = 'M0_Citidilato'
                           elif NM == 'M0_Timidilato':
                              COMPARE = 'M0_Adenilato'
                           elif NM == 'M0_Citidilato':
                              COMPARE = 'M0_Gualinato'
                           print('-------------')
                           PAIRS_DICT[NM] = str("M1_"+ str(COMPARE)[3:])
                           template_ob = bpy.data.objects[COMPARE].select_set(True)
                           bpy.ops.object.duplicate( {COMPARE : template_ob, COMPARE : [template_ob]}, linked=False)
                           obj = bpy.context.selected_objects[-1]
                           obj.location = (OR_OBJ_LC[0]-12,OR_OBJ_LC[1],OR_OBJ_LC[2])
                           obj.name = ("M1_"+ str(COMPARE)[3:])
                           obj.rotation_euler.y = radians(180)
                           bpy.ops.object.select_all(action = 'DESELECT')
                       print('MOLECULES BASIC PAIRS CREATION               - DONE')
                       return()

def OBJ_ROT_CORRRECTION_PA():
                      for obj in bpy.context.scene.objects:
                          if obj.name.startswith("M1"):
                             obj.rotation_mode = 'XYZ'
                             obj.select_set(True)
                             obj.rotation_euler.y = radians(180)
                             obj.rotation_euler.x = radians(180)
                             bpy.ops.object.transform_apply( rotation = True )
                             bpy.ops.object.select_all(action = 'DESELECT')
                      print('OBJ_ROT_CORRRECTION_PA                       - DONE')
                      return()
                      
                      
def TMP_ROT_CORRRECTION_PA():
                      for obj in bpy.context.scene.objects:
                          if obj.name.startswith("M1"):
                             obj.rotation_mode = 'XYZ'
                             obj.select_set(True)
                             #obj.rotation_euler.y = radians(180)
                             obj.rotation_euler.x = radians(180)
                             bpy.ops.object.transform_apply( rotation = True )
                             bpy.ops.object.select_all(action = 'DESELECT')
                      print('OBJ_ROT_CORRRECTION_PA                       - DONE')
                      return()
                      
 
def OBJ_ROT_CORRRECTION_M0():
                      for obj in bpy.context.scene.objects:
                          if obj.name.startswith("M0"):
                             obj.rotation_mode = 'XYZ'
                             obj.select_set(True)
                             obj.rotation_euler.y = radians(180) 
                             bpy.ops.object.transform_apply( rotation = True )
                             bpy.ops.object.select_all(action = 'DESELECT')
                      print('OBJ_ROT_CORRRECTION_M0                       - DONE')
                      return()

def JOIN_PAIRS():
                      bpy.ops.object.select_all(action = 'DESELECT')
                      for NM in PAIRS_DICT:
                            TD = [NM,PAIRS_DICT[NM]]
                            TD
                            temp_NM = "P_"+ str(TD[0])[2:4] +"_" + str(TD[1])[2:4] +"_0"
                            temp_NM
                            for i in TD:
                               obj = bpy.data.objects[i]
                               obj.select_set(True)
                               #bpy.data.objects[i].select_set(True)
                               bpy.context.view_layer.objects.active = obj
                               bpy.ops.object.join()
                            obj.name = temp_NM
                            bpy.ops.object.transform_apply( rotation = True )
                            bpy.data.collections['BASIC_PAIRS'].objects.link(obj)
                            #bpy.data.collections['Collection'].objects.unlink(obj)
                            bpy.ops.object.select_all(action = 'DESELECT')
                      for col in bpy.context.scene.collection.children: 
                          if str(col.name).startswith("DS_") and len(col.objects[:]) > 0:
                             bpy.data.collections[col.name].objects.unlink(col.objects[0])
                             bpy.data.collections.remove(col)  
                      print('JOIN_PAIRS                                   - DONE')
                      return()

def PAIRS_90_ROT():
                    for obj in bpy.context.scene.objects:
                          if obj.name.startswith("P_"):
                             obj.rotation_mode = 'XYZ'
                             obj.select_set(True)
                             obj.rotation_euler.x = radians(90)
                             bpy.ops.object.transform_apply( rotation = True )
                             bpy.ops.object.select_all(action = 'DESELECT')
                    print('PAIRS_90_ROT                                 - DONE')
                    return()
                    
def LOAD_SEQUENCE(): 
                    global DNA_SQN_PATH 
                    with open(DNA_SQN_PATH, 'r') as filehandle:
                         for line in filehandle:
                             currentSQN = line[:-1] 
                    print('LOAD_SEQUENCE                                - DONE')
                    return(currentSQN)
              
def SQNC_3D_GEN(DNA_SQNC):
                    px = 40
                    pz = 0 
                    b_copy = ''
                    for i, v in enumerate(DNA_SQNC):  
                        py = (2.56 * i)
                        if   v == 'a': 
                          b_copy = 'P__A__T_0'
                        elif v == 'g':
                          b_copy = 'P__G__C_0' 
                        elif v == 't':
                          b_copy = 'P__T__A_0' 
                        elif v == 'c':
                          b_copy = 'P__C__G_0'  
                        template_ob = bpy.data.objects[b_copy].select_set(True)
                        bpy.ops.object.duplicate( {b_copy : template_ob, b_copy : [template_ob]}, linked=False)
                        obj = bpy.context.selected_objects[-1]
                        obj.location = (px,py,pz)
                        obj.name = ( str(b_copy)[0:-1]) + "_" + str(i)
                        obj.rotation_euler.y = radians(36*i) #B-DNA 36
                        #obj.rotation_euler.y = radians(33) #A-DNA 31
                        #obj.rotation_euler.y = radians(30) #z-DNA 9/51
                        bpy.ops.object.select_all(action = 'DESELECT')
                    print('SQNC_3D_GEN                                  - DONE')
                    return()
     
              
#for col in bpy.context.scene.collection.ch
#   if str(col.name).startswith("DS_"):    
#     bpy.data.collections.remove(col) 
#------------------------ GLOBAL VARIABLES MANAGEMENT -----------------------------
tt = vars()
s = 'COLOR'
def slicedict(d, s):
    ''' filter variables name'''
    return {k:v for k,v in tt.items() if str(k).startswith(s) }

def FILT_GLOB_VARDICT_to_DF(tt, s):
    ''' filter global variables dictionary and create dataFrame'''
    data_GV=slicedict(tt, s)
    GV_df = pd.DataFrame.from_dict(data_GV, orient='index' )
    return(GV_df)

#----------------------------------------------------------------------------------
SOURCE_PTH   = str('C:\\Users\\ ... FILE_DNA_OK\\')
SOURCE_FL_NM = ''
DNA_SQNC_INPUT_FLDR = str('C:\\ ....\\03_DNA_SEQUENCE_EX\\')
DNA_SQNC_INPUT_FILE = str('DNA_001_125.txt')
DNA_SQN_PATH        = str(DNA_SQNC_INPUT_FLDR + DNA_SQNC_INPUT_FILE )
#----------------------- JSON INPUT
JSON_BASE_DICT_IMPUT = {'A':'DNA_Deoxyadenosine-monophosphate.json','T':'DNA_Thymidine_monophosphate.json','G':'DNA_DeoxyGuanosine_monophosphate.json','C':'DNA_Deoxycytidylate.json'}
COLLECT_DICT         = {'A': 'DS_Adenilato_0', 'T': 'DS_Timidilato_0' ,'G': 'DS_Gualinato_0' ,'C': 'DS_Citidilato_0'}
obj_M0_name          = ['M0_Adenilato', 'M0_Gualinato', 'M0_Timidilato', 'M0_Citidilato']
FUNCTION_LST= ( "Glob_CLLCT_GEN","LOAD_FILE","DF_GENERATOR","MOLECULE_SFARE_GENERATOR","PREPARE_DF_BOUNDS","CREATE_BOUNDS_1","CREATE_LABELS","BASIC_MOLECULE_JOIN","UNLINK_BASIC_MOLECULE_CL","UNLINK_BASIC_MOLECULE_GN","OBJ_MASS_C_MOVE_G","COLUMN_OBJ","CREATE_PAIRS")
PAIRS_DICT = {}
#PAIRS_DICT = {'M0_Adenilato': 'M1_Timidilato', 'M0_Gualinato': 'M1_Citidilato', 'M0_Timidilato': 'M1_Adenilato', 'M0_Citidilato': 'M1_Gualinato'}
#----------- PER ADATTARE COLORI A BLENDER ------------------
sce = bpy.context.scene
ob  = bpy.context.object
display_device = sce.display_settings.display_device
sce.display_settings.display_device = 'None'

#------------------------------------------------------------
GRAY  = bpy.data.materials.new(name="GRAY"  ) #set new material to variable
GRAY.diffuse_color = (0.39,0.38,0.38,0)

try:
    PROS_DF = PROCESS_ROADMAP_PREP(FUNCTION_LST)
except:
    print("PROCESS_ROADMAP_PREP                         - ERROR ")


try:
    Glob_CLLCT_GEN()
except:
    print("COLLECTIONS CREATION                         - ERROR ")

try:
    for i in JSON_BASE_DICT_IMPUT:
            M_TYPE_0 = i
            MOLECULE_NM_L = str(JSON_BASE_DICT_IMPUT[M_TYPE_0])[4:-5]
            print( str(MOLECULE_NM_L) +" file available:   \t"+ str(os.path.isfile(str(SOURCE_PTH + SOURCE_FL_NM ))).upper())
except:
    print("TEST FILE EXISTANCE                               - FAILED")

print("\n  -----------  \n")
try:
    df_ATM.drop()
    df_BND.drop()
except:
    print("NO PREVIOUS DATAFRAME EXISTING or already droped")

X_PLUS_FACTOR = 0
Y_PLUS_FACTOR = 0
M_TYPE_0      = 0

for i in JSON_BASE_DICT_IMPUT:
     index_l = (list(JSON_BASE_DICT_IMPUT.keys()).index(i) )
     M_TYPE_0 = i
     if str(M_TYPE_0)=='A':
        COLLECT_NM = 'DS_Adenilato_0'
     elif str(M_TYPE_0)=='T':
        COLLECT_NM = 'DS_Timidilato_0'
     elif str(M_TYPE_0)=='G':
        COLLECT_NM = 'DS_Gualinato_0'
     elif str(M_TYPE_0)=='C':
        COLLECT_NM = 'DS_Citidilato_0'
     print("PREPARING MOLECULE:  "+COLLECT_NM)
     MOLECULE_NM_L = str(JSON_BASE_DICT_IMPUT[M_TYPE_0])[4:-5]
     #print( MOLECULE_NM_L )
     SOURCE_FL_NM = JSON_BASE_DICT_IMPUT[M_TYPE_0]
     #str(SOURCE_PTH + SOURCE_FL_NM )
     X_PLUS_FACTOR = (int(list(JSON_BASE_DICT_IMPUT.keys()).index(i) )*25)
     try:
         JS_LST_DATA = LOAD_FILE(str(SOURCE_PTH) + SOURCE_FL_NM)
     except:
         print("LOAD PROBLEM                                 - FAILED - Check the input file")
     try:
         df_ATM,df_LEGAMI = DF_GENERATOR(JS_LST_DATA)
     except:
         print("DATAFRAME CREATION PROBLEM - FAILED - Check Json internal structure ")
     try:
         MOLECULE_SFARE_GENERATOR(df_ATM,X_PLUS_FACTOR,COLLECT_NM)
     except:
         print("MOLECULE_SFARE_GENERATOR                     - FAILED")
     try:
         df_BND = PREPARE_DF_BOUNDS(df_ATM,df_LEGAMI)
     except:
         print("BOUNDS PREPARE GENERATION - FASE 1           - FAILED - Check join Datagrame")
     try:
         CREATE_BOUNDS_1(df_BND,X_PLUS_FACTOR,COLLECT_NM)
     except:
         print("BOUNDS GENERATION - FASE 2                   - FAILED")
     bpy.ops.object.select_all(action='DESELECT')
     try:
         CREATE_LABELS(X_PLUS_FACTOR,MOLECULE_NM_L,index_l)
     except:
         print("CREATE_LABELS                                - FAILED")
     print("\n  -----------  \n")

try:
    print("----- BASIC MOLECULE AGGREGATTION -----")
except:
    print(" ")

#----- BASIC MOLECULE AGGREGATTION
CLLCT_LIST = list(COLLECT_DICT.values())

try:
    BASIC_MOLECULE_JOIN(CLLCT_LIST)
except:
    print("MOLECULE AGGREGATION                         - FAILED")

try:
    UNLINK_BASIC_MOLECULE_CL(CLLCT_LIST)
except:
    print("UNLINK_BASIC_MOLECULE_CL                     - FAILED")

try:
    UNLINK_BASIC_MOLECULE_GN(CLLCT_LIST)
except:
    print("UNLINK_BASIC_MOLECULE_GN                     - FAILED")

try:
    OBJ_MASS_C_MOVE_G()
except:
    print("OBJ_MASS_C_MOVE_G                            - FAILED")

try:
    OBJ_ROTATE_G()
except:
    print("OBJ_ROTATE_G                                 - FAILED")

try:
    OBJ_ROT_CORRRECTION()
except:
    print("OBJ_ROT_COR                                  - FAILED")
#------------------
try:
    OBJ_MASS_C_MOVE_G()
except:
    print("OBJ_MASS_C_MOVE_G                            - FAILED")

try:
    COLUMN_OBJ()
except:
    print("COLUMN_OBJ                                   - FAILED")

try:
    CREATE_PAIRS()
except:
    print("CREATE_PAIRS()                               - FAILED")

try:
    OBJ_MASS_C_MOVE_G()
except:
    print("OBJ_MASS_C_MOVE_G                            - FAILED") 
    
try:
    OBJ_ROT_CORRRECTION_PA()
except:
    print("OBJ_ROT_COR                                  - FAILED")

try:
    JOIN_PAIRS()
except:
    print("JOIN_PAIRS                                  - FAILED")

try:
    OBJ_MASS_C_MOVE_G()
except:
    print("OBJ_MASS_C_MOVE_G                            - FAILED")

try:
    PAIRS_90_ROT()
except:
    print("PAIRS_90_ROT                                 - FAILED")
    
try:
    OBJ_MASS_C_MOVE_G()
except:
    print("OBJ_MASS_C_MOVE_G                            - FAILED")

#-------------------------------
try:
    DNA_SQNC = LOAD_SEQUENCE()
except:
    print("LOAD_SEQUENCE                                - FAILED")
    
try: 
    SQNC_3D_GEN(DNA_SQNC)
except:
    print("LOAD_SEQUENCE                                - FAILED")
    
    
    
    
    