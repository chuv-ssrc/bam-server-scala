package scala.controllers.management

import org.scalatestplus.play._
import play.api.libs.json._
import play.api.test.Helpers._
import utils.DbUtils._

import scala.setup.{WithDatabase, WithToken}


/**
  * Test UsersSamplesController.
  */
class UsersSamplesSpec extends PlaySpec with OneAppPerSuite with WithToken with WithDatabase {

  val db = dbContext.db

  "Add an attribution" should {

    "add rows to the users_samples table" in {
      val body = Json.parse(s"""
        { "users_samples": [{"sample": "$testSample3", "username": "$testUser"}, {"sample": "$testSample4", "username": "$testUser"}] }
      """)
      val count0 = countRows(db, "users_samples")
      val response = route(app, FakeAdminRequest(PUT, "/users_samples").withJsonBody(body)).get
      status(response) mustBe OK
      val count = countRows(db, "users_samples")
      count must equal(count0 + 2)
    }

    "ignore adding users_samples that already exist" in {
      val body = Json.parse(s"""
        { "users_samples": [{"sample": "$testSample3", "username": "$testUser"}, {"sample": "$testSample4", "username": "$testUser"}] }
      """)
      val count0 = countRows(db, "users_samples")
      val response = route(app, FakeAdminRequest(PUT, "/users_samples").withJsonBody(body)).get
      status(response) mustBe OK
      val count = countRows(db, "users_samples")
      count must equal(count0)
    }

    "fail if a sample is not known" in {
      val body = Json.parse(s"""
        { "users_samples": [{"sample": "XXXX", "username": "$testUser"}] }
      """)
      val count0 = countRows(db, "users_samples")
      val response = route(app, FakeAdminRequest(PUT, "/users_samples").withJsonBody(body)).get
      status(response) mustBe INTERNAL_SERVER_ERROR
      contentAsString(response) must startWith ("Could not find 'XXXX' in database.")
      val count = countRows(db, "users_samples")
      count must equal(count0)
    }

    "fail if a username is not known" in {
      val body = Json.parse(s"""
        { "users_samples": [{"sample": "$testSample3", "username": "XXXX"}] }
      """)
      val count0 = countRows(db, "users_samples")
      val response = route(app, FakeAdminRequest(PUT, "/users_samples").withJsonBody(body)).get
      status(response) mustBe INTERNAL_SERVER_ERROR
      contentAsString(response) must startWith ("Could not find 'XXXX' in database.")
      val count = countRows(db, "users_samples")
      count must equal(count0)
    }

  }

  "Delete an attribution" should {

    "remove rows from the users_samples table" in {
      val body = Json.parse(s"""
        { "users_samples": [{"sample": "$testSample3", "username": "$testUser"}, {"sample": "$testSample4", "username": "$testUser"}] }
      """)
      val count0 = countRows(db, "users_samples")
      val response = route(app, FakeAdminRequest(DELETE, "/users_samples").withJsonBody(body)).get
      status(response) mustBe OK
      val count = countRows(db, "users_samples")
      count must equal(count0 - 2)
    }

    "ignore deleting non-existent users_samples" in {
      val body = Json.parse(s"""
        { "users_samples": [{"sample": "$testSample3", "username": "$testUser"}, {"sample": "$testSample4", "username": "$testUser"}] }
      """)
      val count0 = countRows(db, "users_samples")
      val response = route(app, FakeAdminRequest(DELETE, "/users_samples").withJsonBody(body)).get
      status(response) mustBe OK
      val count = countRows(db, "users_samples")
      count must equal(count0)
    }

    "fail if a sample is not known" in {
      val body = Json.parse(s"""
        { "users_samples": [{"sample": "XXXX", "username": "$testUser"}] }
      """)
      val count0 = countRows(db, "users_samples")
      val response = route(app, FakeAdminRequest(DELETE, "/users_samples").withJsonBody(body)).get
      status(response) mustBe INTERNAL_SERVER_ERROR
      contentAsString(response) must startWith ("Could not find 'XXXX' in database.")
      val count = countRows(db, "users_samples")
      count must equal(count0)
    }

    "fail if a username is not known" in {
      val body = Json.parse(s"""
        { "users_samples": [{"sample": "otherSample", "username": "XXXX"}] }
      """)
      val count0 = countRows(db, "users_samples")
      val response = route(app, FakeAdminRequest(DELETE, "/users_samples").withJsonBody(body)).get
      status(response) mustBe INTERNAL_SERVER_ERROR
      contentAsString(response) must startWith ("Could not find 'XXXX' in database.")
      val count = countRows(db, "users_samples")
      count must equal(count0)
    }

  }

}
