# CWSI Calculation with Python and Flask API
This repository contains a Python implementation of the Crop Water Stress Index (CWSI) calculation along with a Flask API to expose the functionality as a web service.
This repository contains the Flask application file (main.py), dataset folder (excel_files), backend folder contains all the calculating functions and the jupyter notebook for references (CWSI_calculation.ipynb)

**Table of Contents**
- _Introduction_
- _Installation_
- _Usage_

**Introduction :**
The Crop Water Stress Index (CWSI) is an important indicator used in agricultural research to evaluate water stress in crops. This repository provides a Python library that calculates the CWSI based on preserved data.

**Installation :**
To install the CWSI calculation library and the Flask API, follow these steps:

1. Clone the repository:
```
git clone https://github.com/hohadang99/CWSI-calculation-with-Python-and-FlaskAPI.git
```

2. Create a virtual environment (optional but recommended):
```
python3 -m venv venv
source venv/bin/activate
```

3. Install the required dependencies:
```
pip3 install -r requirements.txt
```

Requirements file including these libraries:
```
numba
numpy==1.21.1
pandas
matplotlib
flask
flask_bootstrap
```

**Usage :**
To calculate the CWSI using the command-line interface, use the following steps:

1. Start Flask app:
```
python app.py
```

2. Go to the Flask address in the command line (E.g. http://127.0.0.1:3000) - The result is shown below:   
![image](https://github.com/hohadang99/CWSI-calculation-with-Python-and-FlaskAPI/assets/40363911/1f016899-80c1-4936-adbc-9395a19e34a1)

3. Click on **Upload file** button and upload data file (the "data.xlsx" file in "excel_files" folder):
![image](https://github.com/hohadang99/CWSI-calculation-with-Python-and-FlaskAPI/assets/40363911/7c57546e-d6b9-45dd-a4aa-f299d6d09ae6)

4. After upload file, choose treatment type, start DOY and end DOY (DOY: Date of the year):
![image](https://github.com/hohadang99/CWSI-calculation-with-Python-and-FlaskAPI/assets/40363911/692f710e-052e-4c2a-a8de-0e032b87ecb9)

5. After step 4, choose the start hour, end hour and hour for CWSI value:
![image](https://github.com/hohadang99/CWSI-calculation-with-Python-and-FlaskAPI/assets/40363911/bfcb0452-c7f0-4136-b406-7261d9b50e12)

6. Showing results:
![image](https://github.com/hohadang99/CWSI-calculation-with-Python-and-FlaskAPI/assets/40363911/6c1f0587-157d-4121-9302-42d77d68be9c)
Click **back to insert value** to start over from step 1
