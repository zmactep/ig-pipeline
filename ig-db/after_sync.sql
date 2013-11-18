#### after syncdb:

drop table django.igtools_train;
drop table django.igtools_predict;
drop table django.igtools_simplecluster;
drop table django.igtools_report;
drop table django.igtools_cutregion;
drop table django.igtools_manifest;
drop table django.igstorage_storageitem;

drop table ig.auth_user_groups;
drop table ig.auth_group_permissions;
drop table ig.auth_group;
drop table ig.auth_user_user_permissions;
drop table ig.auth_permission;
drop table ig.django_admin_log;
drop table ig.django_content_type;
drop table ig.django_session;
drop table ig.django_site;
drop table ig.messages_message;
drop table ig.messages_dialog;
drop table ig.registration_registrationprofile;
drop table ig.auth_user;

insert into ig.igtools_manifest(`tool_name`, `manifest`) values('Train', '{
	"output": [
		{
			"name": "train.libsvm",
			"type": "binary",
			"pipelined": false 
		},
		{
			"name": "train_nominal.arff",
			"type": "text",
			"pipelined": false 
		},
		{
			"name": "model.model",
			"type": "model",
			"pipelined": true 
		}
	]
}'); 

insert into ig.igtools_manifest(`tool_name`, `manifest`) values('Predict', '{
	"output": [
		{
			"name": "debug_prediction.txt",
			"type": "text",
			"pipelined": false 
		},
		{
			"name": "debug_prediction_avg.txt",
			"type": "text",
			"pipelined": false 
		},
		{
			"name": "predict.libsvm",
			"type": "binary",
			"pipelined": false 
		},
		{
			"name": "predict_nominal.arff",
			"type": "text",
			"pipelined": false 
		},
		{
			"name": "predict_nominal_fixed.arff",
			"type": "text",
			"pipelined": false 
		},
		{
			"name": "prediction.txt",
			"type": "text",
			"pipelined": false 
		},
		{
			"name": "read_names.txt",
			"type": "text",
			"pipelined": false 
		},
		{
			"name": "results.kabat",
			"type": "kabat",
			"pipelined": true 
		},
		{
			"name": "results_pic.txt",
			"type": "text",
			"pipelined": false 
		}
	]
}'); 


insert into ig.igtools_manifest(`tool_name`, `manifest`) values('CutRegion', '{
	"output": [
		{
			"name": "output.fasta",
			"type": "fasta",
			"pipelined": true 
		}
	]
}'); 