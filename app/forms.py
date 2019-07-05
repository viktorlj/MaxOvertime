from flask_wtf import FlaskForm
from wtforms import StringField, PasswordField, BooleanField, SubmitField, IntegerField
from wtforms.validators import DataRequired, NumberRange

class NumberOfSamplesForm(FlaskForm):
    numberofsamples = IntegerField('Antal tidpunkter', validators=[DataRequired(), NumberRange(min=2)])
    samplename = StringField('Testnummer', validators=[DataRequired()])
    submit = SubmitField('Fortsätt')

