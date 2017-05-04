package scala.controllers.management

import org.scalatestplus.play._
import play.api.libs.json._
import play.api.test.Helpers._
import utils.DbUtils._

import scala.setup.{WithDatabase, WithToken}


/**
  * Test SamplesController.
  */
class AppsSpec extends PlaySpec with OneAppPerSuite with WithToken with WithDatabase {

  "Add an app" should {

    "add one rows to the app table" in {
      val body = Json.parse("""
        { "iss": "testApp", "keyFile": "/", "description": "none" }
      """)
      val db = dbContext.db
      val count0 = countRows(db, "apps")
      val response = route(app, FakeAdminRequest(PUT, "/apps").withJsonBody(body)).get
      status(response) mustBe OK
      val count = countRows(db, "apps")
      count must equal(count0 + 1)
    }

  }

  "Delete an app" should {

    "remove one rows from the apps table" in {
      val body = Json.parse("""
        { "iss": "testApp" }
      """)
      val db = dbContext.db
      val count0 = countRows(db, "apps")
      val response = route(app, FakeAdminRequest(DELETE, "/apps").withJsonBody(body)).get
      status(response) mustBe OK
      val count = countRows(db, "apps")
      count must equal(count0 - 1)
    }

  }

}
