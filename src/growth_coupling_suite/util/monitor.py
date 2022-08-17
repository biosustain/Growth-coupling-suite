# from weakref import WeakKeyDictionary
from pickle import dumps

class Monitor():
    """
    Monitors objects and checks changes
    """
    def __init__(self):
        # self.objects = WeakKeyDictionary()
        self.objects = {}
        
    def is_changed(self, obj_name, obj):
        current_state = dumps(obj, -1)
        if obj_name in self.objects:
            is_changed = self.objects[obj_name] != current_state
                
        self.objects[obj_name] = current_state
        
        return is_changed
    
    def add_object(self, obj_name, obj):
        self.objects[obj_name] = dumps(obj, -1)
        
        
    
        
        
    