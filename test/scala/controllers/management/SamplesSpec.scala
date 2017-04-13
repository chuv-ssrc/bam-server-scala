package scala.controllers.management

import org.scalatestplus.play._
import play.api.libs.json._
import play.api.test.Helpers._
import utils.DbUtils._

import scala.setup.{WithDatabase, WithToken}


/**
  * Test IndexController.
  */
class SamplesSpec extends PlaySpec with OneAppPerSuite with WithToken with WithDatabase {

  "Add a sample" should {

    "add rows to the samples table" in {
      val body = Json.parse("""
        { "samples": [{"name": "A", "filename": "/"}, {"name": "B", "filename": "/"}] }
      """)
      val db = dbContext.db
      val count0 = countRows(db, "samples")
      val response = route(app, FakeAdminRequest(PUT, "/samples").withJsonBody(body)).get
      status(response) mustBe OK
      val count = countRows(db, "samples")
      count must equal(count0 + 2)
    }

  }

  "Delete a samples" should {

    "remove rows from the samples table" in {
      val body = Json.parse("""
        { "samples": ["A","B"] }
      """)
      val db = dbContext.db
      val count0 = countRows(db, "samples")
      val response = route(app, FakeAdminRequest(DELETE, "/samples").withJsonBody(body)).get
      status(response) mustBe OK
      val count = countRows(db, "samples")
      count must equal(count0 - 2)
    }

  }

}
