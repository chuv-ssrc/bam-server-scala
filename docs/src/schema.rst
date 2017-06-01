

Local database
==============

Here we describe the local bam-server database.
It links together app, user and sample identifers to set individual permissions.


Schema
------

Here is the complete db schema:

.. code:: sql

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


Fields description
------------------
