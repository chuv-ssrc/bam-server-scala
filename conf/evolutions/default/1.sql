# -- Database schema

# --- !Ups

CREATE TABLE `apps` (
    `id` INTEGER(11) NOT NULL PRIMARY KEY AUTO_INCREMENT,
    `iss` VARCHAR(255) NOT NULL,
    `keyFile` VARCHAR(255) NOT NULL,
    `description` VARCHAR(255) DEFAULT NULL,
    `isActive` TINYINT(1) NOT NULL DEFAULT 1
);

CREATE TABLE `users` (
    `id` INTEGER(11) NOT NULL PRIMARY KEY AUTO_INCREMENT,
    `app_id` INTEGER(11) NOT NULL,
    `username` VARCHAR(255) NOT NULL,
    `group` VARCHAR(255) DEFAULT NULL,
    `isActive` TINYINT(1) DEFAULT 0,
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
(1, 'https://jdelafon.eu.auth0.com/', 'auth0.cer', 'Using a public certificate in .cer format', 1),
(2, 'test', 'id_rsa.pub', 'Using only a public key in .pem format', 1)
;

INSERT INTO `users`(`app_id`,`username`,`isActive`,`isAdmin`) VALUES
(1, 'testuser', 1, 0),
(1, 'admin'   , 1, 1),
(1, 'a@test.com', 1, 0),
(1, 'b@test.com', 0, 0),
(1, 'julien.delafontaine@yandex.com', 1, 1)
;

INSERT INTO `samples`(`name`,`filename`,`isActive`) VALUES
('sample1', 'test1.bam', 1),
('sample2', 'test2.bam', 1),
('notherekey', 'nothere.bam', 1),
('inactivekey', 'inactive.bam', 0)
;

INSERT INTO `users_samples`(`user_id`,`sample_id`) VALUES
(2, 1),
(2, 2),
(2, 3),
(2, 4),
(3, 1),
(3, 2),
(3, 3),
(3, 4)
;


# --- !Downs

DROP TABLE `bam`;
