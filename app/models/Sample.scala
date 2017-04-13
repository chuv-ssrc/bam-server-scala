package models

import play.api.libs.json.Json


case class Sample (
  name: String,
  filename: String,
  project: Option[String] = None,
  hash: Option[String] = None,
  description: Option[String] = None,
  isOnDisk: Option[Boolean] = None,
  isActive: Boolean = false
)

object Sample {
  implicit val sampleFormat = Json.format[Sample]
}