create database ig;
use ig;
create table ig_backend_tasks (
	id INT NOT NULL AUTO_INCREMENT,
	params VARCHAR(1024) CHARACTER SET utf8,
	result VARCHAR(1024) CHARACTER SET utf8,
	status ENUM('ok', 'failed', 'in progress'),
	PRIMARY KEY (id)
)
