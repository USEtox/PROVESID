import requests
import logging
import json
import time

class OPSIN:
    def __init__(self):
        self.base_url = "https://opsin.ch.cam.ac.uk/opsin/"
        self.responses = {200: "SUCCESS", 404: "FAILURE", 500: "Internal server error"}

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
            "status": "",
            "message": "",
            "inchi": "",
            "stdinchi": "",
            "stdinchikey": "",
            "smiles": ""
        }
    
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
            if res["status"] == "SUCCESS":
                results.append(res)
            else:
                logging.warning(f"Failed to get ID for {iupac_name}: {res['message']}")
                results.append(res)
            time.sleep(pause_time)
        return results