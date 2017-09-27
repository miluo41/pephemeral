from flask import Flask, render_template, request
from wtforms import Form, FloatField, validators,TextField, widgets, SelectField
from compute import compute
app = Flask(__name__)
from flaskexample import views
