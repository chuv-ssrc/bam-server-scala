# -- Database schema

# --- !Ups

CREATE TABLE `bam` (
    `id` INTEGER(11) NOT NULL PRIMARY KEY AUTO_INCREMENT,
    `filename` VARCHAR(255) NOT NULL,
    `key` VARCHAR(255) NOT NULL,
    `hash` VARCHAR(255) DEFAULT NULL,
    `description` VARCHAR(255) DEFAULT NULL,
    `ondisk` TINYINT(1) DEFAULT NULL,
    `active` TINYINT(1) NOT NULL DEFAULT 0
);

INSERT INTO `bam` VALUES (1, 'test.bam', 'testkey', NULL, NULL, NULL, 1);
INSERT INTO `bam` VALUES (2, 'nothere.bam', 'notherekey', NULL, NULL, NULL, 1);
INSERT INTO `bam` VALUES (3, 'inactive.bam', 'inactivekey', NULL, NULL, NULL, 0);

# --- !Downs

DROP TABLE `bam`;
