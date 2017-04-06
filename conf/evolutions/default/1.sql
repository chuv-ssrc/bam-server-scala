# -- Database schema

# --- !Ups

CREATE TABLE `apps` (
    `id` INTEGER(11) NOT NULL PRIMARY KEY AUTO_INCREMENT,
    `iss` VARCHAR(255) NOT NULL,
    `name` VARCHAR(255) NOT NULL,
    `description` VARCHAR(255) DEFAULT NULL,
);

CREATE TABLE `users` (
    `id` INTEGER(11) NOT NULL PRIMARY KEY AUTO_INCREMENT,
    `app_id` INTEGER(11) DEFAULT NULL,
    `username` VARCHAR(255) NOT NULL,
    `isActive` TINYINT(1) DEFAULT 0,
    `isAdmin` TINYINT(1) DEFAULT 0,
    FOREIGN KEY (`app_id`) REFERENCES `apps`(`id`),
);

CREATE TABLE `bam` (
    `id` INTEGER(11) NOT NULL PRIMARY KEY AUTO_INCREMENT,
    `sample` VARCHAR(255) NOT NULL,
    `filename` VARCHAR(255) NOT NULL,
    `hash` VARCHAR(255) DEFAULT NULL,
    `description` VARCHAR(255) DEFAULT NULL,
    `ondisk` TINYINT(1) DEFAULT NULL,
    `active` TINYINT(1) NOT NULL DEFAULT 0
);

CREATE TABLE `users_bam` (
    `id` INTEGER(11) NOT NULL PRIMARY KEY AUTO_INCREMENT,
    `user_id` INTEGER(11) NOT NULL,
    `bam_id` INTEGER(11) NOT NULL,
    FOREIGN KEY (`user_id`) REFERENCES `users`(`id`),
    FOREIGN KEY (`bam_id`) REFERENCES `bam`(`id`),
);


INSERT INTO `apps` VALUES (1, 'https://jdelafon.eu.auth0.com/', 'auth0', NULL);

INSERT INTO `users` VALUES (1,1, 'testuser', 1, 0);
INSERT INTO `users` VALUES (2,1, 'admin', 1, 1);
INSERT INTO `users` VALUES (3,1, 'julien.delafontaine@yandex.com', 1, 0);

INSERT INTO `bam` VALUES (1, 'testkey', 'test.bam', NULL, NULL, NULL, 1);
INSERT INTO `bam` VALUES (2, 'notherekey', 'nothere.bam', NULL, NULL, NULL, 1);
INSERT INTO `bam` VALUES (3, 'inactivekey', 'inactive.bam', NULL, NULL, NULL, 0);

INSERT INTO `users_bam` VALUES (1, 1, 1);
INSERT INTO `users_bam` VALUES (2, 1, 2);
INSERT INTO `users_bam` VALUES (3, 1, 3);
INSERT INTO `users_bam` VALUES (4, 2, 1);
INSERT INTO `users_bam` VALUES (5, 2, 2);
INSERT INTO `users_bam` VALUES (6, 2, 3);
INSERT INTO `users_bam` VALUES (7, 3, 1);
INSERT INTO `users_bam` VALUES (8, 3, 2);
INSERT INTO `users_bam` VALUES (9, 3, 3);


# --- !Downs

DROP TABLE `bam`;
