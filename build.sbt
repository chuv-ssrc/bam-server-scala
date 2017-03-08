import NativePackagerHelper._

name := """bam-server"""
version := "1.0-SNAPSHOT"
scalaVersion := "2.11.7"

lazy val root = (project in file(".")).enablePlugins(PlayScala)

libraryDependencies ++= Seq(
  jdbc,
  cache,
  ws,
  "org.scalatestplus.play" %% "scalatestplus-play" % "1.5.0" % Test,
  "com.github.samtools" % "htsjdk" % "2.1.1",
  "mysql" % "mysql-connector-java" % "5.1.36",
  "org.mockito" % "mockito-all" % "1.9.5"
)

resolvers += "scalaz-bintray" at "http://dl.bintray.com/scalaz/releases"

libraryDependencies += filters

// Add non-default directories to build (/target/universal/) using NativePackageHelper
mappings in Universal ++= directory("scripts")