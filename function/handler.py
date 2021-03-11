import json
import os
from flask import Request

from rdkit.Chem import Descriptors, AllChem as Chem, DataStructs
import numpy as np
import joblib
import pickle
from .selected_targets import o_dict, pn_dict, th_dict

function_root = os.environ.get("function_root")

N_BITS = 1024
INPUT_DIR = "/home/app/chembl_mcp_models"


# load models and scalers ----------------------------------------------------

models = {}
scalers = {}
for target in pn_dict.keys():
    # load models
    model_path = f"{INPUT_DIR}/models/{target}/{target}_conformal_prediction_model"
    models[target] = joblib.load(model_path)
    # load scalers
    scaler_path = f"{INPUT_DIR}/scalers/{target}_scaler.pkl"
    scalers[target] = pickle.load(open(scaler_path, "rb"))

# ----------------------------------------------------------------------------


def pred_category(p0, p1, significance):
    if (p0 >= significance) & (p1 >= significance):
        return "both"
    if (p0 >= significance) & (p1 < significance):
        return "inactive"
    if (p0 < significance) & (p1 >= significance):
        return "active"
    else:
        return "empty"


def predict(descriptors):
    predictions = []
    for target in models:
        # load scaler
        scaler = scalers[target]
        # predict
        X = np.column_stack(
            (
                scaler.transform(np.array(descriptors[:6]).reshape(1, -1)),
                descriptors[-1].reshape(1, -1),
            )
        )
        pred = models[target].predict(X)
        # get p values
        p0 = float(pred[:, 0])
        p1 = float(pred[:, 1])
        # format output for a single prediction
        res = {
            "target_chemblid": target,
            "organism": o_dict[target],
            "pref_name": pn_dict[target],
            "70%": pred_category(p0, p1, 0.3),
            "80%": pred_category(p0, p1, 0.2),
            "90%": pred_category(p0, p1, 0.1),
            "threshold": th_dict[target]
        }
        predictions.append(res)
    return predictions


def calc_descriptors(rdmol):
    fp = Chem.GetMorganFingerprintAsBitVect(
        rdmol, radius=2, nBits=N_BITS, useFeatures=False
    )
    np_fp = np.zeros(N_BITS)
    ecfp = DataStructs.ConvertToNumpyArray(fp, np_fp)
    logp = Descriptors.MolLogP(rdmol)
    mwt = Descriptors.MolWt(rdmol)
    rtb = Descriptors.NumRotatableBonds(rdmol)
    hbd = Descriptors.NumHDonors(rdmol)
    hba = Descriptors.NumHAcceptors(rdmol)
    tpsa = Descriptors.TPSA(rdmol)
    return [logp, mwt, rtb, hbd, hba, tpsa, np_fp]


def handle(req: Request):
    """handle a request to the function.

    Your response is immediately passed to the caller, unmodified.
    This allows you full control of the response, e.g. you can set
    the status code by returning a tuple (str, int). A detailed
    description of how responses are handled is found here:

    http://flask.pocoo.org/docs/1.0/quickstart/#about-responses

    Args:
        req (Request): Flask request object
    """
    predictions = []
    # load molecule from smiles and calculate fp
    mol = Chem.MolFromSmiles(req.get_json()["smiles"])
    if mol:
        # calc descriptors
        descs = calc_descriptors(mol)
        # predict
        predictions = predict(descs)

    return json.dumps(predictions)
