# CWSI Calculation with Python and Flask API
This repository contains a Python implementation of the Crop Water Stress Index (CWSI) calculation along with a Flask API to expose the functionality as a web service.
This repository contains the Flask application file (main.py), the dataset folder (excel_files), the backend folder contains all the calculating functions, and the jupyter notebook for references (CWSI_calculation.ipynb)

**Table of Contents**
- _Introduction_
- _Installation_
- _Usage_

**Introduction:**
The Crop Water Stress Index (CWSI) is an important indicator used in agricultural research to evaluate water stress in crops. This repository provides a Python library that calculates the CWSI based on preserved data.

**Installation:**  
To install the CWSI calculation library and the Flask API, follow these steps:

1. Clone the repository:
```
git clone https://github.com/lhomanh/CWSI-Calculation-with-Python-and-FlaskAPI.git
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

**Usage:**
To calculate the CWSI using the command-line interface, use the following steps:

1. Start Flask app:
```
python app.py
```

2. Go to the Flask address in the command line (E.g. http://127.0.0.1:3000) - The result is shown below:   
![image](https://github.com/lhomanh/CWSI-Calculation-with-Python-and-FlaskAPI/assets/60242356/64a0e411-ee68-47d9-8cf9-2c2374465660)


3. Click on the **Upload file** button and upload the data file (the "data.xlsx" file in the "excel_files" folder):
![image](https://github.com/lhomanh/CWSI-Calculation-with-Python-and-FlaskAPI/assets/60242356/9a485f61-ccf4-45f1-8004-aa17b8598fe5)


4. After uploading the file, choose treatment type, start DOY, and end DOY (DOY: Date of the year):
![image](https://github.com/lhomanh/CWSI-Calculation-with-Python-and-FlaskAPI/assets/60242356/7d9b0717-b5d2-4566-a89f-181913e74fe6)


5. After step 4, choose the start hour, end hour, and hour for CWSI value:
![image](https://github.com/lhomanh/CWSI-Calculation-with-Python-and-FlaskAPI/assets/60242356/adf0373a-2797-421c-aafa-30afa4461331)


6. Showing results:
![image](https://github.com/lhomanh/CWSI-Calculation-with-Python-and-FlaskAPI/assets/60242356/a1f8e6f9-c515-4681-9df2-4adcf9821a2d)

Click **back to insert value** to start over from step 1
