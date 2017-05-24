# -- Database schema

# --- !Ups

CREATE TABLE `apps` (
    `id` INTEGER(11) NOT NULL PRIMARY KEY AUTO_INCREMENT,
    `iss` VARCHAR(255) NOT NULL,
    `key` TEXT NOT NULL,
    `algorithm` VARCHAR(255) DEFAULT 'RS256',
    `description` VARCHAR(255) DEFAULT NULL,
    `isActive` TINYINT(1) NOT NULL DEFAULT 1
);

CREATE TABLE `users` (
    `id` INTEGER(11) NOT NULL PRIMARY KEY AUTO_INCREMENT,
    `app_id` INTEGER(11) NOT NULL,
    `username` VARCHAR(255) NOT NULL,
    `group` VARCHAR(255) DEFAULT NULL,
    `isActive` TINYINT(1) DEFAULT 1,
    `isAdmin` TINYINT(1) DEFAULT 0,
    FOREIGN KEY (`app_id`) REFERENCES `apps`(`id`)
);

CREATE TABLE `samples` (
    `id` INTEGER(11) NOT NULL PRIMARY KEY AUTO_INCREMENT,
    `name` VARCHAR(255) NOT NULL,
    `filename` VARCHAR(255) NOT NULL,
    `project` VARCHAR(255) DEFAULT NULL,
    `hash` VARCHAR(255) DEFAULT NULL,
    `description` VARCHAR(255) DEFAULT NULL,
    `isOndisk` TINYINT(1) DEFAULT NULL,
    `isActive` TINYINT(1) NOT NULL DEFAULT 1
);

CREATE TABLE `users_samples` (
    `id` INTEGER(11) NOT NULL PRIMARY KEY AUTO_INCREMENT,
    `user_id` INTEGER(11) NOT NULL,
    `sample_id` INTEGER(11) NOT NULL,
    `isActive` TINYINT(1) NOT NULL DEFAULT 1,
    FOREIGN KEY (`user_id`) REFERENCES `users`(`id`),
    FOREIGN KEY (`sample_id`) REFERENCES `samples`(`id`)
);


INSERT INTO `apps` VALUES
(1, 'https://jdelafon.eu.auth0.com/',
'-----BEGIN RSA PUBLIC KEY-----\nMIIBIjANBgkqhkiG9w0BAQEFAAOCAQ8AMIIBCgKCAQEAtWWKxPv9vsWdRR/hcmJF\nsQjjUrMs/OsVstyNJXwmWuhl3lNIZwwEDoJbnE9IKPyizyNwbnB9FmJnClCboUeP\nbkuIrDM63+S+PtX/SQ9YI5yDxz+88dRYT86WP23wcWMO3txV2GAu62RVGSl48ZJP\nSyu94NBIiZOO5oDJpWDInhZphiMQ3u/rEwlVxVMt0CTTInfl4iX0sCtymD2y6M38\nVrQwHOzSddFrbI58t4Rfal4SttwdmXONRnj7mrgl5G6v7IHEa/HOrlT1rSLOMBKz\nOfmZy+bdlt5zrx3Adfzgn1BC6DGlG3Y9QYMOPpXjbzRO3rv9Fl5bRJyn5Ih82Cey\ndQIDAQAB\n-----END RSA PUBLIC KEY-----',
'RS256', 'Using a public certificate in .cer format', 1),
(2, 'testapp', '-----BEGIN RSA PUBLIC KEY-----\nMFwwDQYJKoZIhvcNAQEBBQADSwAwSAJBAK7ttYaE/1ldsb0OJQDQhhDWqwuFWIyt\nxgYIJH1HYA4UpA/Nm24fERIA1xi2Pomep6VTnQ/ThFP5hn2NyITwCIsCAwEAAQ==\n-----END RSA PUBLIC KEY-----',
'RS256', 'Using a public key in .pem format', 1),
(3, 'testapp2', 'secretHMACkey',
'HS256', 'Using a shared secret key', 1)
;

INSERT INTO `users`(`id`,`app_id`,`username`,`isActive`,`isAdmin`) VALUES
(1, 1, 'admin@test.com', 1, 1),
(2, 1, 'test@test.com', 1, 0),
(3, 2, 'test@test.com', 1, 0),
(4, 3, 'test@test.com', 1, 0)
;

INSERT INTO `samples`(`id`,`name`,`filename`,`isActive`) VALUES
(1, 'sample1', 'test1.bam', 1),
(2, 'sample2', 'test2.bam', 1),
(3, 'sample3', 'nothere.bam', 1),
(4, 'sample4', 'nothere.bam', 1),
(5, 'notherekey', 'nothere.bam', 1),
(6, 'inactivekey', 'inactive.bam', 0)
;

INSERT INTO `users_samples`(`user_id`,`sample_id`) VALUES
(1, 1),(1, 2), (1, 5),(1, 6),
(2, 1),        (2, 5),(2, 6)
;


# --- !Downs

DROP TABLE `bam`;
