# -- Database schema

# --- !Ups

CREATE TABLE `bam` (
    `id` INTEGER(11) NOT NULL PRIMARY KEY AUTO_INCREMENT,
    `username` VARCHAR(255) NOT NULL,
    `sample` VARCHAR(255) NOT NULL,
    `filename` VARCHAR(255) NOT NULL,
    `hash` VARCHAR(255) DEFAULT NULL,
    `description` VARCHAR(255) DEFAULT NULL,
    `ondisk` TINYINT(1) DEFAULT NULL,
    `active` TINYINT(1) NOT NULL DEFAULT 0
);

INSERT INTO `bam` VALUES (1, 'testuser', 'testkey', 'test.bam', NULL, NULL, NULL, 1);
INSERT INTO `bam` VALUES (2, 'testuser', 'notherekey', 'nothere.bam', NULL, NULL, NULL, 1);
INSERT INTO `bam` VALUES (3, 'testuser', 'inactivekey', 'inactive.bam', NULL, NULL, NULL, 0);

# --- !Downs

DROP TABLE `bam`;
