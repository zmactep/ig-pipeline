from django.db import models
import datetime
from django.utils import timezone

# Create your models here.

class Poll(models.Model):
	question = models.CharField(max_length=200)
	pub_date = models.DateTimeField('date published')

	def __str__(self):
		return self.question

	def was_published_recently(self):
		return self.pub_date >= timezone.now() - datetime.timedelta(days=1)

class Choise(models.Model):
	poll = models.ForeignKey(Poll)
	choise_text = models.CharField(max_length=200)
	votes = models.IntegerField(default=0)

	def __str__(self):
		return self.choise_text
