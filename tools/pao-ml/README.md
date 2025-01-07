# PAO-ML training

These are the scripts for training an equivariant PAO-ML model. To create a new model, follow these
steps:

1. Create and activate a [virtual Python environment](https://docs.python.org/3/tutorial/venv.html):

   ```
   python3 -m venv ./venv
   source ./venv/bin/activate
   ```

1. Install the required Python packages:

   ```
   pip3 install torch e3nn scipy
   ```

1. Train a new model:

   ```
   ./pao-train.py --kind=H training_data1.pao training_data2.pao...
   ```

1. Retrain the model on some more taining data:

   ```
   ./pao-retrain.py --model="DZVP-MOLOPT-GTH-PAO4-H.pt" training_data100.pao training_data101.pao...
   ```

1. Validate the model against test data:

   ```
   ./pao-validate.py --model="DZVP-MOLOPT-GTH-PAO4-H.pt" test_data1.pao test_data2.pao...
   ```

1. Use the model in a CP2K run by setting the `PAO_MODEL_FILE` keyword in the kind section:

   ```
   &KIND H
     BASIS_SET DZVP-MOLOPT-GTH
     PAO_BASIS_SIZE 4
     PAO_MODEL_FILE DZVP-MOLOPT-GTH-PAO4-H.pt
   &END KIND
   ```
