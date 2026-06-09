# Helper: writes a minimal connection.yaml in tempdir, invokes
# load_connection() against it, then cleans up. Keeps every test
# isolated from the developer's real config.
load_connection_test <- function(deployment, allow_optin = FALSE) {

    tmp_dir <- withr::local_tempdir()
    config_path <- file.path(tmp_dir, "connection.yaml")
    yaml::write_yaml(
        list(
            default_deployment = "dev3",
            credentials_source = "env"
        ),
        config_path
    )

    load_connection(
        deployment = deployment,
        config_path = config_path,
        allow_optin = allow_optin
    )
}


test_that("panel_deployment falls back to the global default", {

    cfg <- list(default_deployment = "dev3", overrides = list())

    expect_equal(
        panel_deployment("any-figure", config = cfg),
        "dev3"
    )
    expect_equal(
        panel_deployment("any-figure", facet = "any-facet", config = cfg),
        "dev3"
    )
})

test_that("built-in registry routes the FR-030 panels to dev4", {

    cfg <- list(default_deployment = "dev3", overrides = list())

    expect_equal(
        panel_deployment("fig01-overview", "structures", config = cfg),
        "dev4"
    )
    expect_equal(
        panel_deployment("fig01-overview", "panel_e", config = cfg),
        "dev4"
    )
    expect_equal(
        panel_deployment("fig01-overview", "ramp_conflict", config = cfg),
        "dev4"
    )
    expect_equal(
        panel_deployment("tab02-ramp-comparison", config = cfg),
        "dev4"
    )
    expect_equal(
        panel_deployment("fig01-overview", config = cfg),
        "dev3"
    )
})

test_that("user-config overrides win over the built-in registry", {

    cfg <- list(
        default_deployment = "dev3",
        overrides = list(
            `fig01-overview` = list(structures = "dev5")
        )
    )

    expect_equal(
        panel_deployment("fig01-overview", "structures", config = cfg),
        "dev5"
    )
    # The override only matches its specific facet — fall back otherwise.
    expect_equal(
        panel_deployment("fig01-overview", "panel_e", config = cfg),
        "dev4"
    )
})

test_that("user-config panel-level default wins over the global default", {

    cfg <- list(
        default_deployment = "dev3",
        overrides = list(
            `fig07-supplement` = list(default = "dev4")
        )
    )

    expect_equal(
        panel_deployment("fig07-supplement", config = cfg),
        "dev4"
    )
    # A facet-specific lookup that doesn't match anything still falls
    # back to the panel-level default before the global default.
    expect_equal(
        panel_deployment("fig07-supplement", "unknown", config = cfg),
        "dev4"
    )
})

test_that("load_connection refuses dev5 (reserved for integrated build)", {

    expect_error(
        load_connection_test("dev5"),
        regexp = "reserved"
    )
})

test_that("load_connection gates prod / dev2 behind allow_optin", {

    expect_error(
        load_connection_test("prod"),
        regexp = "excluded"
    )
    expect_error(
        load_connection_test("dev2"),
        regexp = "excluded"
    )

    record_prod <- load_connection_test("prod", allow_optin = TRUE)
    expect_equal(record_prod$name, "prod")
    expect_equal(record_prod$role, "opt-in")
})

test_that("load_connection resolves dev3 and dev4 records", {

    record3 <- load_connection_test("dev3")
    expect_equal(record3$name, "dev3")
    expect_equal(record3$db_port, 5403L)
    expect_equal(record3$role, "figures-default")

    record4 <- load_connection_test("dev4")
    expect_equal(record4$name, "dev4")
    expect_equal(record4$db_port, 5404L)
    expect_equal(record4$role, "figures-structures")
})

test_that("PGPORT_<DEPLOYMENT> env override wins over the table port", {

    withr::with_envvar(
        c(PGPORT_DEV3 = "55403"),
        {
            record <- load_connection_test("dev3")
            expect_equal(record$db_port, 55403L)
        }
    )
})
