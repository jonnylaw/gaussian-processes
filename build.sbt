lazy val commonSettings = Seq(
  scalaVersion := "2.12.7",
  name := "gaussian-processes",
  organization := "com.github.jonnylaw",
  version := "0.1.0",
  scalacOptions ++= Seq(
    "-encoding", "UTF-8",   // source files are in UTF-8
    "-deprecation",         // warn about use of deprecated APIs
    "-unchecked",           // warn about unchecked type parameters
    "-feature",             // warn about misused language features
    "-language:higherKinds",// allow higher kinded types without `import scala.language.higherKinds`
    "-Xlint",               // enable handy linter warnings
//    "-Xfatal-warnings",     // turn compiler warnings into errors
    "-Ypartial-unification", // allow the compiler to unify type constructors of different arities
    "-language:implicitConversions" // allow implicit conversion of DLM -> DGLM
  ),
  credentials += Credentials(Path.userHome / ".sbt" / ".credentials"),
  licenses := Seq("APL2" -> url("http://www.apache.org/licenses/LICENSE-2.0.txt")),
  homepage := Some(url("https://jonnylaw.github.io/gaussian-processes")),
  scmInfo := Some(
    ScmInfo(
      url("https://github.com/jonnylaw/gaussian-processes"),
      "scm:git@github.com:jonnylaw/gaussian-processes.git"
    )
  ),
  developers := List(
    Developer(id="1", name="Jonny Law", email="law.jonny@googlemail.com", url=url("https://jonnylaw.github.io/blog"))
  ),
  publishMavenStyle := true,
  publishTo := {
    val nexus = "https://oss.sonatype.org/"
    if (isSnapshot.value)
      Some("snapshots" at nexus + "content/repositories/snapshots")
    else
      Some("releases"  at nexus + "service/local/staging/deploy/maven2")
  }
)

//ensimeScalaVersion in ThisBuild := "2.12.4"

lazy val micrositeSettings = Seq(
  micrositeName := "Gaussian Processes",
  micrositeDescription := "Gaussian Processes",
  micrositeBaseUrl := "/gaussian-processes",
  micrositeDocumentationUrl := "/gaussian-processes/docs",
  micrositeGithubOwner := "jonnylaw",
  micrositeGithubRepo := "gaussian-processes",
  micrositeImgDirectory := (resourceDirectory in Compile).value / "figures",
  micrositeCssDirectory := (resourceDirectory in Compile).value / "styles",
  micrositeHighlightTheme := "solarized-dark",
  micrositeGithubToken := sys.env.get("GITHUB_TOKEN"),
  micrositeShareOnSocial := false,
  micrositePushSiteWith := GitHub4s,
  micrositeCDNDirectives := microsites.CdnDirectives(
    jsList = List(
      "https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.5/MathJax.js?config=TeX-MML-AM_CHTML"
    ))
)

lazy val core = (project in file("core"))
  .settings(
    commonSettings,
    resolvers += Resolver.bintrayRepo("cibotech", "public"),
    libraryDependencies ++= Seq(
      "org.scalanlp"             %% "breeze"     % "0.13.2",
      "com.github.fommil.netlib" % "all"         % "1.1.2",
      "com.cibo"                 %% "evilplot"   % "0.6.3",
      "com.47deg"                %% "github4s"   % "0.19.0",
      "org.typelevel"            %% "cats-core"  % "1.5.0",
      "org.scalatest"            %% "scalatest"  % "3.0.5"  % "test",
      "org.scalacheck"           %% "scalacheck" % "1.13.4" % "test"
    ),
  )

lazy val docs = (project in file("docs"))
  .settings(moduleName := "docs")
  .settings(micrositeSettings: _*)
  .enablePlugins(MicrositesPlugin)
  .enablePlugins(TutPlugin)
  .dependsOn(core)

lazy val examples = (project in file("examples"))
  .settings(
    libraryDependencies ++= Seq(
      "com.nrinaudo"           %% "kantan.csv"         % "0.4.0",
      "com.nrinaudo"           %% "kantan.csv-generic" % "0.4.0",
      "com.github.jonnylaw"    %% "bayesian_dlms"      % "0.4.0",
      "com.github.nscala-time" %% "nscala-time"        % "2.20.0"
    )
  )
  .dependsOn(core)
