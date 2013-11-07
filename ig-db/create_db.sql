create database ig;
create database django;
use ig;
create table ig_backend_tasks (
	id INT NOT NULL AUTO_INCREMENT,
	params VARCHAR(4096) CHARACTER SET utf8,
	result VARCHAR(4096) CHARACTER SET utf8,
	status ENUM('ok', 'failed', 'in progress'),
	PRIMARY KEY (id)
);


#### after syncdb:

drop table django.igsnooper_taskrequest;
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