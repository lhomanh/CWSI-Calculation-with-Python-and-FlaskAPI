from backend.func_new import *
from lib2to3.pgen2.pgen import DFAState
from flask import Flask, request, redirect, url_for, render_template, flash,Blueprint
from flask_bootstrap import Bootstrap
from werkzeug.utils import secure_filename


main = Flask(__name__)
Bootstrap(main)
main.secret_key = b'_5#y2L"F4Q8z\n\xec]/'
OUTPUT_DIR = 'static'   
main.config['UPLOAD_FOLDER'] = OUTPUT_DIR
upload_dir =  './excel_files/'


@main.route('/')
def index():
    return render_template('/base/index.html')

@main.route('/return_upload')
def return_upload():
    return render_template('/base/upload.html')

@main.route('/upload')
def upload():
    return render_template('/base/upload.html')

@main.route('/profile')
def profile():
    return render_template('/base/profile.html')


@main.route('/upload',methods=['POST'])
def upload_file():
    if request.method == 'POST':
        if 'file' not in request.files:
            print('No file part')
            return redirect(request.url)
        file = request.files['file']
        if file.filename == '':
            print('No selected file')
            return redirect(request.url)
        if file:
            global filename_
            filename_ = file.filename
        return render_template('/base/upload_success.html')


@main.route('/upload_success',methods=['GET','POST'])
def upload_success():     
    if request.method == 'POST':
        full_name = upload_dir+filename_  
        global treatment
        treatment = request.form['treatment']
        global start_DOY
        start_DOY = request.form['start_DOY']
        start_DOY = int(start_DOY)
        global end_DOY
        end_DOY = request.form['end_DOY']
        end_DOY = int(end_DOY)
        print("NEXT")
        return render_template("/base/upload_success_2nd.html")
        
    return render_template('/base/upload_success.html')



@main.route('/upload_success_2nd',methods=['GET','POST'])
def upload_success_2nd():     
    if request.method == 'POST':
        full_name = upload_dir+filename_  
        global treatment
        global start_DOY
        global end_DOY
        start_hour = request.form['start_hour']
        start_hour = int(start_hour)
        end_hour = request.form['end_hour']
        end_hour = int(end_hour)
        select_hour = request.form['select_hour']
        select_hour = int(select_hour)
        print("NEXT")
        plot_buttons(full_name, start_DOY, end_DOY, treatment, start_hour, end_hour, select_hour)
        return render_template("/base/result.html")
        
    return render_template('/base/upload_success.html')


if __name__ == "__main__":
    main.run(host="0.0.0.0", port=3000, threaded=False,debug=True)
