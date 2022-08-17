# =============================================================================
# Module for using the KEGG API 
# =============================================================================
import requests
import re
import json
from os.path import dirname, abspath, isdir
from os import listdir, makedirs
from io import StringIO
import pandas as pd

# get directory of module
dir_module = dirname(abspath(__file__))


# %%


class KeggAPI():
    
    def __init__(self):
        # initialize KEGG API
        
        
        # %% globals
        self.dir = dir_module
        self.kegg_query_address = "http://rest.kegg.jp/find/{0}/{1}/{2}" # for querying a search word
        self.kegg_get_address = "http://rest.kegg.jp/get/{0}/{1}" # for retrieving a given database entry
        self.kegg_link_address = "http://rest.kegg.jp/link/{0}/{1}" # for linking an entry to a database
        self.kegg_list_address = "http://rest.kegg.jp/list/{0}" # for listing descriptions of an entry
        
    
        # save queries in this session
        self.library_filename = "kegg_query_library"       
        # load library       
        self.query_library = {}
        self.query_dir_format = self.dir + "/query_library/kegg/{0}/" # 0: database
        self.query_file_format = "{0}_q_kegg_{1}.json" # 0: query text; 1: database
        self.query_file_format_option = "{0}_{2}_q_kegg_{1}.json" # 0: query text; 1: database; 2: option
        self.load_query_library()
        

    def query_kegg(self, query, database, query_address, option='', requery=False):
        # request query from KEGG API
        
        # check option
        if not(option):
            option = ""
            query_key = query
        else:
            query_key = query + "_" + option
        
        
        # check if database is present
        if database not in self.query_library:
            # create database library
            self.create_database_library(database)
            
        # check if query result already exists
        if not(requery) and (query_key in self.query_library[database]):
            # load result
            r_dict = self.query_library[database][query_key]
            text = r_dict["text"]
            status = r_dict["status_code"]
            
        else:  
            # query kegg
            # r = requests.get(self.kegg_query_address.format(database, query, option))
            r = requests.get(query_address)

            # save query in library
            # self.query_library[database][query] = r
        
            
            if r.status_code == 200:
                # save response 
                self.query_library[database][query_key] = {"query_text": query,
                                                       "text": r.text,
                                                       "status_code": r.status_code,
                                                       "option": option}
                # extract response 
                text = r.text
                status = r.status_code
                
            else:
                # save response 
                self.query_library[database][query_key] = {"query_text": query,
                                                       "text": None,
                                                       "status_code": r.status_code,
                                                       "option": option}
                # extract response 
                text = None
                status = r.status_code
            
            # save query
            self.save_query(query_key, database, self.query_library[database][query_key])
            
        return text, status
    
    
    def get(self, database, query, option='', prioritize=False, best_guess=False):
        # database: KEGG database to be queried
                        # compound
                        
                        
        # query: query text
        # option: options text
        # prioritize: determine and return only the closest match 
        
           
        # query kegg
        query_address = self.kegg_query_address.format(database, query, option)
        query_text, status_code = self.query_kegg(query, database, query_address, option=option)
        
        #  query "compound" database
        if database=="compound":
            # evaluate text from KEGG API query
            compounds = self.evaluate_query_text(query_text, database)
            
            # check if compound was found
            if not(compounds["ID"]):
                # no compounds found on KEGG
                if best_guess:
                    # make a guess for the compound name, gradually erase first letters
                    number_erase_letters = 7 # erase maximal n letters
                    for i in range(number_erase_letters):
                        query_truncated = query[(i+1):]
                        # query kegg
                        query_address = self.kegg_query_address.format(database, query_truncated, "")
                        query_text, status_code = self.query_kegg(query_truncated, database, query_address)
                        # evaluate text
                        compounds = self.evaluate_query_text(query_text, database)
                        # compound found?
                        if not(compounds["ID"]):
                            pass
                        else:
                            query = query_truncated # further use truncated query word
                            break
                    
            if prioritize and len(compounds["ID"])>0:
                
                # check all synonyms and rate match grade for all
                
                
                query_letters = "".join(re.findall("[a-z]|[A-Z]", query))
                # determine the exact match to the query word, disregard case
                for i in range(len(compounds["synonyms"])):
                    if query.lower() in (synonym.lower() for synonym in compounds["synonyms"][i]):
                        compound_match = {}
                        compound_match["synonyms"] = compounds["synonyms"][i]
                        compound_match["ID"] = compounds["ID"][i]
                        # return match
                        compound_match["kegg_status_code"] = status_code
                        return compound_match
                    else:
                        match_grade = []
                 
                        for synonym in compounds["synonyms"][i]:
                            synonym_letters = "".join(re.findall("[a-z]|[A-Z]", synonym))
                            if query_letters.lower() in synonym_letters.lower():
                                match_grade.append(len(query_letters) / len(synonym_letters))
                            else:
                                match_grade.append(0)
                                
                        compounds["match"].append(max(match_grade))
                                
                    
                # if no match was found get the closest match
                # match_list = []
                # match_grades = compounds["match"]
                # for item in compounds["match"]: match_list.append(item)
                # # in case there is no match return all synonyms
                # if not match_list:
                #     return compounds
                
                if any(compounds["match"]):
                    best_match = compounds["match"].index(max(compounds["match"]))
                    compound_match = {}
                    compound_match["synonyms"] = compounds["synonyms"][best_match]
                    compound_match["ID"] = compounds["ID"][best_match]
                    compound_match["match"] = compounds["match"][best_match]
                    # return closest match
                    compound_match["kegg_status_code"] = status_code
                    return compound_match
                else:
                    pass
    
            
            # return results
            compounds["kegg_status_code"] = status_code
            return compounds
     
        
    def get_entry(self, entry_key, option=''):
        """retrieve specific database entry"""
        
        
        database = 'entries'
        
        
        # query kegg
        query_address = self.kegg_get_address.format(entry_key, option)
        text, status = self.query_kegg(entry_key, database, query_address)

            
        return text, status
    
    
    
    def link_entry(self, entry_key, database, option='', requery=False) -> pd.DataFrame:
        """retrieve specific database entry"""
        
        # query kegg
        query_address = self.kegg_link_address.format(database, entry_key)
        text, status = self.query_kegg(entry_key, database, query_address, requery=requery)
   
            
        # convert text to dataframe   
        column_names = ['Input', 'Output']
        if text:
            # entries = [e.split[for e in text.split('\n')
            s = StringIO(text)
            return_frame = pd.read_csv(s, sep='\t', header=None, names=column_names)
        else:
            return_frame = pd.DataFrame(columns=column_names)
            
            
        return return_frame, status        
       
    
    
    def list_entry(self, entry_key, option='', requery=False) -> list:
        """list descriptions of a kegg entry"""
        
        database = 'list'

        # query kegg
        query_address = self.kegg_list_address.format(entry_key)
        text, status = self.query_kegg(entry_key, database, query_address, requery=requery)

        # convert text to dataframe  
        column_names = ['Input', 'Output']
        if text:
            # entries = [e.split[for e in text.split('\n')
            s = StringIO(text)
            return_frame = pd.read_csv(s, sep='\t', header=None, names=column_names)
            
        else:
            return_frame = pd.DataFrame(columns=column_names)
            
        return return_frame, status
        
        
    # %%
        
    def evaluate_query_text(self, text, database):
        # evaluate text returned by KEGG API
        
        # compound database
        if database == "compound":
            compounds = {}
            compounds["ID"] = []
            compounds["synonyms"] = []
            compounds["match"] = []
            # process compound queries
            # identifiy compounds
            spl = re.split("\t|\n",text)
            for i in range(len(spl)):
                # define content            
                if spl[i].find("cpd:") != -1:   
                    # define a new compound           
                    cmp = spl[i][4:len(spl[i])]
                    # get and save all synonyms         
                    # save in arrays
                    compounds["synonyms"].append(spl[i+1].split("; "))
                    compounds["ID"].append(cmp)
                    
            return compounds
        
        
    def save_query_library(self):
        # save query library to the hard drive
        with open(self.dir + "/" + self.library_filename + ".json", "w") as f:
            json.dump(self.query_library, f)
            
        
    def save_query(self, query_text, database, response):
        # correct query text
        query_text_corrected = query_text.replace(':', '_')
        # save a single query
        try:
            with open(self.query_dir_format.format(database) \
                      + self.query_file_format.format(query_text_corrected, database), "w") as f:
                json.dump(response, f)
        except:
            print('KEGG query result could not be saved.')
        
     
    def create_database_library(self, database):
        """Create folder and dict for new database library"""
        
        dirname = self.query_dir_format.format(database)
        if not isdir(dirname):
            makedirs(dirname)
        else:
            print('Library', database, 'already exists')
            
        if database not in self.query_library:
            self.query_library[database] = {}
        else:
            print('Library', database, 'already in memory')
 
        
        
    def load_query_library(self, override=False):
        # load query library
        
        # default library design
        standard_databases = ['list', "entries"]
        # query_library_default = {"compound": {}
        #                          }
                                 
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
        
        query_databases = [d
                            for d in listdir(self.query_dir_format.format(''))
                            if isdir(self.query_dir_format.format(d))
                            ]        
        
        for database in query_databases+standard_databases:
            self.query_library[database] = {}
            # get all files
            dirname = self.query_dir_format.format(database)
            if not isdir(dirname):
                makedirs(dirname)
                continue
            filenames = listdir(dirname)
            for filename in filenames:
                if self.query_file_format.format("", database) in filename:
                    # load query
                    with open(dirname + filename, "r") as f:
                        query_dict = json.load(f)
                        self.query_library[database][query_dict["query_text"]] \
                            = query_dict
            
        
                                   
        # if isfile(filename):
        #     # load existing libary
        #     with open(filename, "r") as f:
        #         self.query_library = json.load(f)
            
        #     # check if all database entries exist
        #     keys_exist = list(self.query_library.keys())
        #     for key in query_library_default.keys():
        #         if not(key in keys_exist):
        #             self.query_library[key] = {}
                               
        # else:
        #     self.query_library = query_library_default
                
        
                                         
        
        
        
