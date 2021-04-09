from flask_wtf import FlaskForm
from wtforms import StringField, SubmitField, IntegerField, FileField, MultipleFileField
from wtforms.validators import DataRequired, NumberRange


class RegistrationForm(FlaskForm):
    freq= IntegerField('Frequency', validators=[DataRequired(), NumberRange(min=0, max=99)])
    submit = SubmitField('Submit')


class FileForm(FlaskForm):
    fileinput = MultipleFileField('Input', validators=[DataRequired()])
            

