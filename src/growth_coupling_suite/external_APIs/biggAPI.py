# BIGG API for querying BIGG database

import requests
import json
from time import perf_counter as pc
from time import sleep
from os.path import dirname, abspath, isdir
from os import mkdir, makedirs, listdir

# get directory of module
dir_module = dirname(abspath(__file__))

# 
class BiggAPI():
    def __init__(self):
        # get directory of module and query library
        self.dir = dir_module
        
        # initialize timer
        self._last_query_time = pc()
        # minimum time between two BIGG calls [sec]
        self.min_time_diff = 0.1
        # BIGG query address
        self.bigg_query_address = "http://bigg.ucsd.edu/api/v2/{0}/{1}/{2}"
        
        
        # save queries in this session
        self.library_filename = "bigg_query_library"
        
        # load library and setup library directory                 
        self.query_library = {}
        self.query_dir_format = self.dir + "/query_library/bigg/{0}/{1}/" # 0: database; 1: model
        self.query_file_format = "{0}_q_bigg_{1}_{2}.json" # 0: query text; 1: database; 2: model
        self.load_query_library(self.dir + "/" + self.library_filename + ".json")
        
        
        
        
        
    def get(self, query_text, database, model="universal"):

        # query BIGG Database and return result
        if database == "metabolites":
            return self._get_metabolite(query_text, model)
        
        else:
            raise TypeError("Unknown database type: " + database)
        
       
        

    
    
    def _get_metabolite(self, query_text, model):
        # metabolite specific query address
        database_name = "metabolites"
        
        # check if query result already exists
        if model in self.query_library[database_name]:
            if query_text in list(self.query_library[database_name][model].keys()):
                # load result
                return self.query_library[database_name][model][query_text]
        else:
            self.query_library[database_name][model] = {}
        
        res = self._query_bigg(self.bigg_query_address.format(model,
                                                          database_name,
                                                          query_text))

        
        if self._check_query(res):
            # query successful
            res_return = res.json()
            
        else:
            # query failed
            res_return = {}   
         
        # save query text
        res_return["query_text"] = query_text
        # save results
        self.query_library[database_name][model][query_text] = res_return
        self.save_query(query_text, database_name, model, res_return)
            
        return res_return
         
    
    
    def _query_bigg(self, query_address):
        # check timing
        
        # check time of last query
        time_diff = pc()-self._last_query_time
        
        if time_diff < self.min_time_diff:
            # wait
            sleep(self.min_time_diff-time_diff)
              
        # query bigg database
        res = requests.get(query_address)
         
        # reset timer
        self._last_query_time =  pc()
         
        return res
       
    
    def save_query_library(self):
        # save query library to the hard drive
        with open(self.dir + "/" + self.library_filename + ".json", "w") as f:
            json.dump(self.query_library, f)
     
        
    def save_query(self, query_text, database, model, response):
        # save a single query
        dirname = self.query_dir_format.format(database, model)
        if not(isdir(dirname)):
            # create model directory
            mkdir(dirname)
        
        with open(dirname \
                  + self.query_file_format.format(query_text, database, model), "w") as f:
            json.dump(response, f) 
        
        
    def load_query_library(self, override=False):
        # load query library
        
    
        
        # default library design
        query_databases = ["metabolites"]

                                 
        # check if library is already loaded
        if not(self.query_library):
            # library is empty, load library from file if exists
            pass
        else:
            print("Warning: Query library already in memory!")
            if override:
                print("Load new library nevertheless...")
            else:
                print("Abort")
                return
               
        # create base folders if not existent
        if not isdir(self.query_dir_format.format('')):
            makedirs(self.query_dir_format.format(''))
        
        for database in query_databases:
            # init database
            self.query_library[database] = {}
            # get directory for database
            dirname = self.query_dir_format.format(database, "")[:-1]
            # check if path exists
            if not isdir(dirname):
                makedirs(dirname)
                continue
            # get model specific queries
            models = [folder for folder in listdir(dirname) if isdir(dirname + folder)]
            for model in models:            
                self.query_library[database][model] = {}
                # get all files
                dirname = self.query_dir_format.format(database, model)
                filenames = listdir(dirname)
                for filename in filenames:
                    if self.query_file_format.format("", database, model) in filename:
                        # load query
                        with open(dirname + filename, "r") as f:
                            query_dict = json.load(f)
                            self.query_library[database][model][query_dict["query_text"]] \
                                = query_dict
        

    
    def _check_query(self, res):
        # check query result
        check_success = True
        if res.status_code != 200 or not(res.ok):
            # query failed
            check_success = False
        
        return check_success