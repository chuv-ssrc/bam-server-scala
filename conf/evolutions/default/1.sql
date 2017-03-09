# Database schema

# --- !Ups

CREATE TABLE `bam` (
    `id` INTEGER(11) NOT NULL AUTO_INCREMENT,
    `filename` VARCHAR(255) NOT NULL,
    `key` VARCHAR(255) NOT NULL,
    `hash` VARCHAR(255) DEFAULT NULL,
    `description` VARCHAR(255) DEFAULT NULL,
    `ondisk` VARCHAR(255) DEFAULT NULL,
    `active` TINYINT(1) DEFAULT 1,
    PRIMARY KEY (id)
);

# --- !Downs

DROP TABLE `bam`;
