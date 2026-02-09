import requests
import logging
import json
import time
from .cache import cached, clear_opsin_cache, get_opsin_cache_info
from py2opsin import py2opsin

class OPSIN:
    def __init__(self, use_cache: bool = True):
        self.base_url = "https://opsin.ch.cam.ac.uk/opsin/"
        self.responses = {200: "SUCCESS", 404: "FAILURE", 500: "Internal server error"}
        self.use_cache = use_cache

    def clear_cache(self):
        """Clear the cache for all OPSIN methods"""
        clear_opsin_cache()
    
    def get_cache_info(self):
        """Get information about the current cache state"""
        return get_opsin_cache_info()

    @cached(service='opsin')
    def get_id(self, iupac_name: str, timeout=30):
        """
        Returns the SMILES for a given IUPAC name. The code is adapted from IUPAC WorlFair book:
        https://iupac.github.io/WFChemCookbook/tools/opsin_api_jn.html
        """
        apiurl = self.base_url + iupac_name + '.json'
        res = self._empty_res()
        reqdata = requests.get(apiurl, timeout=timeout)
        if reqdata.status_code != list(self.responses.keys())[0]:
            res["status"] = self.responses[reqdata.status_code]
            return res 
        jsondata = reqdata.json()
        res["iupac_name"] = iupac_name
        res["status"] = self.responses[reqdata.status_code]
        res["smiles"] = jsondata["smiles"]
        res["stdinchi"] = jsondata["stdinchi"]
        res["stdinchikey"] = jsondata["stdinchikey"]
        res["inchi"] = jsondata["inchi"]
        return res

    @staticmethod
    def _empty_res():
        """
        create an empty response dictionary of the following format:
        {
            "status": "SUCCESS",
            "message": "",
            "inchi": "InChI=1/C2H2Cl4/c3-1(4)2(5)6/h1-2H",
            "stdinchi": "InChI=1S/C2H2Cl4/c3-1(4)2(5)6/h1-2H",
            "stdinchikey": "QPFMBZIOSGYJDE-UHFFFAOYSA-N",
            "smiles": "ClC(C(Cl)Cl)Cl"
        }
        """
        return {
            "iupac_name": "",
            "status": "",
            "message": "",
            "inchi": "",
            "stdinchi": "",
            "stdinchikey": "",
            "smiles": ""
        }
    
    @cached(service='opsin')
    def get_id_from_list(self, iupac_names: list, timeout=30, pause_time=0.5):
        """
        Returns a list of dictionaries with the SMILES for a given list of IUPAC names.
        """
        results = []
        for iupac_name in iupac_names:
            try:
                res = self.get_id(iupac_name, timeout)
            except requests.RequestException as e:
                logging.error(f"Request failed for {iupac_name}: {e}")
                res = self._empty_res()
                res["status"] = "FAILURE"
                res["iupac_name"] = iupac_name
            if res["status"] == "SUCCESS":
                results.append(res)
            else:
                logging.warning(f"Failed to get ID for {iupac_name}: {res['message']}")
                results.append(res)
            time.sleep(pause_time)
        return results
    
class PYOPSIN:
    """
    This class uses py2opsin (https://github.com/JacksonBurns/py2opsin) package 
    that can be installed via pip:
    pip install py2opsin
    It provides all functionalities as the OPSIN class above -and some more- but works offline 
    and faster. You need to have Java installed on your system.
    Note:
    output_format (str, optional) – One of “SMILES”, “ExtendedSMILES”, 
    “CML”, “InChI”, “StdInChI”, or “StdInChIKey”. Defaults to “SMILES”.
    """
    def __init__(self, jar_fpath = "default"):
        self.jar_fpath = jar_fpath

    def get_smiles(self, iupac_name: str):
        smiles = py2opsin(iupac_name, output_format="SMILES", jar_fpath=self.jar_fpath)
        return smiles
    
    def get_extended_smiles(self, iupac_name: str):
        ext_smiles = py2opsin(iupac_name, output_format="ExtendedSMILES", jar_fpath=self.jar_fpath)
        return ext_smiles
    
    def get_inchi(self, iupac_name: str):
        inchi = py2opsin(iupac_name, output_format="InChI", jar_fpath=self.jar_fpath)
        return inchi
    
    def get_std_inchi(self, iupac_name: str):
        std_inchi = py2opsin(iupac_name, output_format="StdInChI", jar_fpath=self.jar_fpath)
        return std_inchi
    
    def get_std_inchikey(self, iupac_name: str):
        std_inchikey = py2opsin(iupac_name, output_format="StdInChIKey", jar_fpath=self.jar_fpath)
        return std_inchikey
    
    def get_CML(self, iupac_name: str):
        cml = py2opsin(iupac_name, output_format="CML", jar_fpath=self.jar_fpath)
        return cml

    def get_id(self, iupac_name: str):
        res = self._empty_res()
        res["iupac_name"] = iupac_name
        res["smiles"] = self.get_smiles(iupac_name)
        res["extended_smiles"] = self.get_extended_smiles(iupac_name)
        res["inchi"] = self.get_inchi(iupac_name)
        res["stdinchi"] = self.get_std_inchi(iupac_name)
        res["stdinchikey"] = self.get_std_inchikey(iupac_name)
        res["cml"] = self.get_CML(iupac_name)
        if res["smiles"] == "":
            res["status"] = "FAILURE"
        else:
            res["status"] = "SUCCESS"                
        return res

    def get_id_from_list(self, iupac_names: list):
        """
        This function does not need a loop since the py2opsin package handles lists internally.
        """
        res = self.get_id(iupac_names)
        # convert to list of dicts
        results = []
        for i, name in enumerate(iupac_names):
            single_res = self._empty_res()
            single_res["iupac_name"] = name
            single_res["smiles"] = res["smiles"][i]
            single_res["extended_smiles"] = res["extended_smiles"][i]
            single_res["inchi"] = res["inchi"][i]
            single_res["stdinchi"] = res["stdinchi"][i]
            single_res["stdinchikey"] = res["stdinchikey"][i]
            single_res["cml"] = res["cml"][i]
            single_res["status"] = res["status"][i]
            results.append(single_res)
        return results

    @staticmethod
    def _empty_res():
        """
        “SMILES”, “ExtendedSMILES”, 
    “CML”, “InChI”, “StdInChI”, or “StdInChIKey”. Defaults to “SMILES”.
        create an empty response dictionary of the following format:
        {
            "status": "SUCCESS",
            "message": "",
            "inchi": "InChI=1/C2H2Cl4/c3-1(4)2(5)6/h1-2H",
            "stdinchi": "InChI=1S/C2H2Cl4/c3-1(4)2(5)6/h1-2H",
            "stdinchikey": "QPFMBZIOSGYJDE-UHFFFAOYSA-N",
            "smiles": "ClC(C(Cl)Cl)Cl",
            "extended_smiles": "",
            "cml": ""
        }
        """
        return {
            "status": "",
            "message": "",
            "inchi": "",
            "stdinchi": "",
            "stdinchikey": "",
            "smiles": "",
            "extended_smiles": "",
            "cml": ""
        }

