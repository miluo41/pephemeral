from flaskexample import app
from flaskexample import compute
from flask import Flask, render_template, request
from wtforms import Form, FloatField, validators,TextField, widgets, SelectField


class InputForm(Form):
	def value_check(form,field):
		for i in field.data:
			if i not in 'aArRnNdDcCqQeEgGhHiIlLkKmMfFpPsStTwWyYvVxX[];, ':
				raise validators.ValidationError('Unrecognized amino acid abbreviation')
	r = TextField(validators = [validators.InputRequired(message='Please enter your peptide sequence'), \
		validators.Length(min=8,message=('Please enter a longer peptide')),value_check])
	N_ter = SelectField("N-terminal Modification",choices=[('Free','Free'),('Acetylation','Acetylation'),('Glycosylation','Glycosylation'),('Hydroxylation','Hydroxylation')])
	C_ter = SelectField("C-terminal Modification",choices=[('Free','Free'),('Amidation','Amidation'),('Pegylation','Pegylation'),('Propylamidation','Propylamindation')])
	CL = SelectField('Cyclic or Linear',choices=[('Linear','Linear'),('Cyclic','Cyclic')])
	Assay = SelectField('Assay Method',choices=[('in vivo','in vivo'),('in vitro','in vitro')])

@app.route('/',methods=['GET','POST'])
def index():
	form =InputForm(request.form)
	if request.method == 'POST' and form.validate():
		r = form.r.data
		cat_dict={'in_vivo_in_vitro':form.Assay.data,'Linear_cyclic':form.CL.data,
		'N_ter_mod':form.N_ter.data,'C_ter_mod':form.C_ter.data}
		s = compute(r,cat_dict)
		return render_template('view_output.html',form = form, s=s)
	else:
		return render_template('view_input.html',form = form)
