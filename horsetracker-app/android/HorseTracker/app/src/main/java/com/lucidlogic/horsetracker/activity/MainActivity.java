package com.lucidlogic.horsetracker.activity;

import android.os.Bundle;
import android.support.v4.app.Fragment;
import android.view.View;
import android.support.design.widget.NavigationView;
import android.support.v4.view.GravityCompat;
import android.support.v4.widget.DrawerLayout;
import android.support.v7.app.ActionBarDrawerToggle;
import android.support.v7.app.AppCompatActivity;
import android.support.v7.widget.Toolbar;
import android.view.Menu;
import android.view.MenuItem;

import com.lucidlogic.horsetracker.R;
import com.lucidlogic.horsetracker.config.Constants;
import com.lucidlogic.horsetracker.nav.AppNavigationManager;
import com.lucidlogic.horsetracker.view.impl.RacecardFragment;
import com.lucidlogic.horsetracker.view.impl.ResultFragment;
import com.lucidlogic.horsetracker.view.impl.StableFragment;

import timber.log.Timber;

public class MainActivity extends AppCompatActivity
        implements NavigationView.OnNavigationItemSelectedListener {

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        AppNavigationManager.initialize(this);
        setContentView(R.layout.activity_main);
        Toolbar toolbar = (Toolbar) findViewById(R.id.toolbar);
        setSupportActionBar(toolbar);

        DrawerLayout drawer = (DrawerLayout) findViewById(R.id.drawer_layout);
        ActionBarDrawerToggle toggle = new ActionBarDrawerToggle(
                this, drawer, toolbar, R.string.navigation_drawer_open, R.string.navigation_drawer_close);
        drawer.addDrawerListener(toggle);
        toggle.syncState();

        NavigationView navigationView = (NavigationView) findViewById(R.id.nav_view);
        navigationView.setNavigationItemSelectedListener(this);

        addFragment(Constants.RACECARD);

    }

    private void addFragment(String label){
        Fragment frag = AppNavigationManager.findFragment(label);
        if(frag  == null){
            Timber.i(" new instance %s", label);
            switch (label) {
                case Constants.RACECARD: {
                    frag = RacecardFragment.newInstance();
                    break;
                }
                case Constants.RESULT: {
                    frag = ResultFragment.newInstance();
                    break;
                }
                case Constants.STABLE: {
                    frag = StableFragment.newInstance("", "");
                    break;
                }
            }
        }
        AppNavigationManager.addFragment(frag, label);
    }

    @Override
    public void onBackPressed() {
        DrawerLayout drawer = (DrawerLayout) findViewById(R.id.drawer_layout);
        if (drawer.isDrawerOpen(GravityCompat.START)) {
            drawer.closeDrawer(GravityCompat.START);
        } else {
            onBackPressedFragment();
        }
    }

    public void onBackPressedFragment() {

        int count = getFragmentManager().getBackStackEntryCount();

        if (count == 0) {
            super.onBackPressed();
        } else {
            getFragmentManager().popBackStack();
        }
    }

    @Override
    public boolean onCreateOptionsMenu(Menu menu) {
        // Inflate the menu; this adds items to the action bar if it is present.
        getMenuInflater().inflate(R.menu.main, menu);
        return true;
    }

    @Override
    public boolean onOptionsItemSelected(MenuItem item) {
        // Handle action bar item clicks here. The action bar will
        // automatically handle clicks on the Home/Up button, so long
        // as you specify a parent activity in AndroidManifest.xml.
        int id = item.getItemId();

        //noinspection SimplifiableIfStatement
        if (id == R.id.action_settings) {
            Timber.i("action_settings");
            return true;
        } else if(id == R.id.action_other){
            Timber.i("action_other");
            return true;
        }

        return super.onOptionsItemSelected(item);
    }

    @SuppressWarnings("StatementWithEmptyBody")
    @Override
    public boolean onNavigationItemSelected(MenuItem item) {
        // Handle navigation view item clicks here.
        int id = item.getItemId();

        if (id == R.id.nav_racecard) {
            Timber.i("Racecard");
            addFragment(Constants.RACECARD);
        } else if (id == R.id.nav_result) {
            addFragment(Constants.RESULT);
        } else if (id == R.id.nav_stable) {
            Timber.i("My Stable");
            addFragment(Constants.STABLE);
        }

        DrawerLayout drawer = (DrawerLayout) findViewById(R.id.drawer_layout);
        drawer.closeDrawer(GravityCompat.START);
        return true;
    }

    public void forceCrash(View view) {
        throw new RuntimeException("This is a crash");
    }

}
