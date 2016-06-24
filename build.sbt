name := """bam-server"""

version := "1.0-SNAPSHOT"

lazy val root = (project in file(".")).enablePlugins(PlayScala)

scalaVersion := "2.11.7"

libraryDependencies ++= Seq(
  jdbc,
  cache,
  ws,
  "org.scalatestplus.play" %% "scalatestplus-play" % "1.5.0-RC1" % Test,
  "com.github.samtools" % "htsjdk" % "2.1.1",
  "mysql" % "mysql-connector-java" % "5.1.36"

)

resolvers += "scalaz-bintray" at "http://dl.bintray.com/scalaz/releases"
