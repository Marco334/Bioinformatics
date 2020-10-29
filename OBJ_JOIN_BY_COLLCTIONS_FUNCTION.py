import os
import bpy

def OBJ_JOIN_BY_COLLCTIONS(CLLCT_LIST):
                       for COLL_NM in CLLCT_LIST: #Collection list) 
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
                       print('COLLECTION OBJECT JOIN - DONE')                           
                       return()
