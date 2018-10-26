package com.lucidlogic.horsetracker.config;

import android.content.Context;

import com.crashlytics.android.Crashlytics;
import com.lucidlogic.horsetracker.BuildConfig;
import com.lucidlogic.horsetracker.R;
import com.lucidlogic.horsetracker.logging.ReleaseTree;
import com.lucidlogic.horsetracker.model.dto.EntityDTO;
import com.lucidlogic.horsetracker.model.dto.EntityDetailFutureDTO;
import com.lucidlogic.horsetracker.model.dto.EntityDetailHistoricalDTO;
import com.lucidlogic.horsetracker.model.dto.MeetingDTO;
import com.lucidlogic.horsetracker.model.dto.RaceDTO;
import com.lucidlogic.horsetracker.model.dto.RacecardDTO;
import com.lucidlogic.horsetracker.model.dto.RunnerDTO;
import com.lucidlogic.horsetracker.model.dto.UserFollowsEntityDTO;
import com.parse.Parse;
import com.parse.ParseObject;
import com.parse.http.ParseHttpRequest;
import com.parse.http.ParseHttpResponse;
import com.parse.http.ParseNetworkInterceptor;

import java.io.IOException;

import io.fabric.sdk.android.Fabric;
import timber.log.Timber;

/**
 * Created by ian on 06/07/16.
 */
public class AppConfig {

    public static void init(Context context){

        //Parse Init
        ParseObject.registerSubclass(RaceDTO.class);
        ParseObject.registerSubclass(MeetingDTO.class);
        ParseObject.registerSubclass(RacecardDTO.class);
        ParseObject.registerSubclass(RunnerDTO.class);
        ParseObject.registerSubclass(EntityDTO.class);
        ParseObject.registerSubclass(EntityDetailFutureDTO.class);
        ParseObject.registerSubclass(EntityDetailHistoricalDTO.class);
        ParseObject.registerSubclass(UserFollowsEntityDTO.class);

        String appId = context.getString(R.string.parse_app_id);
        String serverUrl = context.getString(R.string.parse_server_url);
        String clientId = context.getString(R.string.parse_client_id);

        Parse.initialize(
            new Parse.Configuration.Builder(context)
                .applicationId(appId)
                .enableLocalDataStore()
                .addNetworkInterceptor(new ParseNetworkInterceptor() {
                    @Override
                    public ParseHttpResponse intercept(Chain chain) throws IOException {
                        ParseHttpRequest req = chain.getRequest();
                        Timber.i(req.getUrl());
                        ParseHttpResponse resp = chain.proceed(req);
                        Timber.i(String.format("resp status code %d", resp.getStatusCode()));
                        return resp;
                    }
                })
                .server(serverUrl)
                .clientKey(clientId)
                .build()
        );

        if (BuildConfig.DEBUG) {
            Timber.plant(new Timber.DebugTree());
        } else {
            Timber.plant(new ReleaseTree());
            //Crashlytics
            Fabric.with(context, new Crashlytics());
        }

    }
}