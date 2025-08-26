import os
import json
import requests
import logging
CASCommonChem_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')

class CASCommonChem:
    """
    calling the CAS Common Chemistry API to get the information for a given CAS RN
    """
    def __init__(self, swagger_file_name='commonchemistry-swagger.json'):
        self.data_folder = CASCommonChem_path
        self.swagger_file_path = os.path.join(self.data_folder, swagger_file_name)
        # load the json file to a dictionary
        with open(self.swagger_file_path, 'r') as f:
            self.swagger = json.load(f)
        host = self.swagger["host"]
        schemes = self.swagger["schemes"][0]
        base_path = self.swagger["basePath"]
        self.base_url = f"{schemes}://{host}{base_path}"
        self.query_url = ["/detail", "/export", "/search"]
        self.responses = {200: "Success", 400: "Invalid Request", 404: "Invalid Request", 500: "Internal Server Error"}
    
    def cas_to_detail(self, cas_rn: str, timeout=30):
        """
        Returns a dictionary with the data for a given CAS RN, the cas number must be a string with or without hyphens.
        A dictionary is returned with the following keys:
        status: the status of the request
        canonicalSmile: the canonical smile
        experimentalProperties: a list of experimental properties
        hasMolfile: a boolean value
        images: a list of images
        inchi: the inchi
        inchiKey: the inchi key
        molecularFormula: the molecular formula
        molecularMass: the molecular mass
        name: the name of the compound
        propertyCitations: a list of property citations
        replacedRns: a list of replaced CAS RN
        rn: the CAS RN
        smile: the smile
        synonyms: a list of synonyms
        uri: the uri
        """
        url = self.base_url + self.query_url[0] + "?cas_rn=" + cas_rn
        res = self._empty_res()
        try:
            response = requests.get(url, timeout=timeout)
            res["status"] = self.responses[response.status_code]
        except:
            res["status"] = "Error"
            return res
        
        if response.status_code != 200:
            return res
        data = response.json()
        for key in data.keys():
            res[key] = data[key]
        return res
    
    def name_to_detail(self, name: str, timeout=30):
        """
        Returns a dictionary with the data for a given name. It works with SMILES too.
        """
        res = self._empty_res()
        url = self.base_url + self.query_url[2] + "?q=" + name
        try:
            response = requests.get(url, timeout=timeout)
        except:
            res["status"] = "Error"
            return res
        # get the CASRN from the response
        res_call = response.json()
        if "count" not in res_call or res_call["count"]==0:
            res["status"] = "Not found"
            return res 
        if res_call["count"]>1:
            logging.warning(f"Multiple compounds found for {name}")
        cas_rn = res_call["results"][0]["rn"]
        return self.cas_to_detail(cas_rn)
    
    def smiles_to_detail(self, smiles: str, timeout=30):
        return self.name_to_detail(smiles, timeout)
    
    @staticmethod
    def _empty_res():
        return {
                "cas_rn": "",
                "status": "",
                "canonicalSmile": "",
                "experimentalProperties": [],
                "hasMolfile": False,
                "images": [],
                "inchi": "",
                "inchiKey": "",
                "molecularFormula": "",
                "molecularMass": "",
                "name": "",
                "propertyCitations": [],
                "replacedRns": [],
                "rn": "",
                "smile": "",
                "synonyms": [],
                "uri": ""
            }
        
        
        
