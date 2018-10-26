package com.lucidlogic.horsetracker.service;

import com.lucidlogic.horsetracker.BuildConfig;
import com.lucidlogic.horsetracker.MockApp;

import org.junit.runner.RunWith;
import org.robolectric.RobolectricGradleTestRunner;
import org.robolectric.annotation.Config;

import static org.mockito.Mockito.mock;

/**
 * Created by ian on 10/07/16.
 */
@RunWith(RobolectricGradleTestRunner.class)
@Config(constants = BuildConfig.class, sdk = 21, application = MockApp.class)
public class ParseTest {

}
