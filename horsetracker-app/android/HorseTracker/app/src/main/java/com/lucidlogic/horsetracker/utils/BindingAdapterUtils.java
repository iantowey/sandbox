package com.lucidlogic.horsetracker.utils;

import android.databinding.BindingAdapter;
import android.support.v4.view.PagerAdapter;
import android.support.v4.view.ViewPager;
import android.widget.Button;

import timber.log.Timber;

/**
 * Created by ian on 28/07/16.
 */
public class BindingAdapterUtils {

    @BindingAdapter("textConditional")
    public static void setOnPageChangeListener(Button btn, boolean val) {
        Timber.i("textConditional");
        btn.setText(val ? "Unfollow" : "Follow");
    }

    @BindingAdapter("setPageChangeListener")
    public static void setOnPageChangeListener(ViewPager viewPager, ViewPager.OnPageChangeListener listener) {
        Timber.i("setOnPageChangeListener");
        viewPager.addOnPageChangeListener(listener);
    }

    @BindingAdapter("setPagerAdapter")
    public static void setOnPageChangeListener(ViewPager viewPager, PagerAdapter pagerAdapter) {
        Timber.i("setPagerAdapter");
        viewPager.setAdapter(pagerAdapter);
    }
}
