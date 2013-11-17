create database ig;
create database django;
create table ig.igbackend_tasks (
	id INT NOT NULL AUTO_INCREMENT,
	params VARCHAR(4096) CHARACTER SET utf8,
	result VARCHAR(4096) CHARACTER SET utf8,
	status ENUM('ok', 'failed', 'in progress'),
	PRIMARY KEY (id)
);
