package scala.controllers.management

import org.scalatestplus.play._
import play.api.libs.json._
import play.api.test.Helpers._
import utils.DbUtils._

import scala.setup.{WithDatabase, WithToken}


/**
  * Test SamplesController.
  */
class SamplesSpec extends PlaySpec with OneAppPerSuite with WithToken with WithDatabase {

  val db = dbContext.db

  "Add a sample" should {

    "add rows to the samples table" in {
      val body = Json.parse("""
        { "samples": [{"name": "A", "filename": "/"}, {"name": "B", "filename": "/"}] }
      """)
      val count0 = countRows(db, "samples")
      val response = route(app, FakeAdminRequest(PUT, "/samples").withJsonBody(body)).get
      status(response) mustBe OK
      val count = countRows(db, "samples")
      count must equal(count0 + 2)
    }

    "fail and insert nothing if one of the samples already exists" in {
      val body = Json.parse("""
        { "samples": [{"name": "A", "filename": "/"}, {"name": "B", "filename": "/"}] }
      """)
      val count0 = countRows(db, "samples")
      val response = route(app, FakeAdminRequest(PUT, "/samples").withJsonBody(body)).get
      status(response) mustBe INTERNAL_SERVER_ERROR
      val count = countRows(db, "samples")
      count must equal(count0)
    }

  }

  "Delete a samples" should {

    "remove rows from the samples table" in {
      val body = Json.parse("""
        { "samples": ["A","B"] }
      """)
      val count0 = countRows(db, "samples")
      val response = route(app, FakeAdminRequest(DELETE, "/samples").withJsonBody(body)).get
      status(response) mustBe OK
      val count = countRows(db, "samples")
      count must equal(count0 - 2)
    }

    "fail and delete nothing if one of the samples does not exist" in {
      val body = Json.parse("""
        { "samples": ["A","B"] }
      """)
      val count0 = countRows(db, "samples")
      val response = route(app, FakeAdminRequest(DELETE, "/samples").withJsonBody(body)).get
      status(response) mustBe INTERNAL_SERVER_ERROR
      val count = countRows(db, "samples")
      count must equal(count0)
    }

  }

}
